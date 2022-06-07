#include "../include/Simulation.h"
#include "../include/ODE.h"
#include "../include/Util.h"
#include "../include/OrbitalDynamics.h"
#include "../include/OrbitalElements.h"

#include "../lib/Eigen/Dense"

#define MAX_HALO_ITERATIONS 10

const double Section::SECTION_DAYS = 1;

Section::Section(EarthMoonSun* f){
    this->dynamics = f;
    this->ode = ODE_RK4<6>();
    this->ode.set_dynamics(f);
    this->ode.recording.set_record_interval(120);
    this->ode.set_timestep(10);
}

void Section::compute_states() {
    this->ode.recording.clear();
    this->ode.set_time(this->t_start);
    this->ode.run(this->initial_state,this->t_final);
    this->times = this->ode.recording.get_times();
    this->states = this->ode.recording.get_states();
    this->times.push_back(ode.get_time());
    this->states.push_back(ode.get_state());
}

double Section::target(const std::array<double,6>& xp){
	this->compute_states();
    this->compute_STM();
    
    const std::array<double,6>& xf = this->states.back();
    
    Matrix<3,4> L;
	Vector<3> b;
    for(uint_fast8_t i = 0; i < 3; i++) {
        for(uint_fast8_t j = 0; j < 3; j++) {
            L[i][j] = this->STM[i][j+3];
        }
        L[i][3] = xf[i + 3];
        b.data[i] = xp[i] - xf[i];
    }

	double delta  = sqrt(b.data[0]*b.data[0] + b.data[1]*b.data[1] + b.data[2]*b.data[2]);

	double stepSize = 1/(1 + delta*5e-4);

    Matrix<4,3> LT = L.get_transpose();
    Matrix<3,3> A = L*LT;
    Matrix<3,3>::solve(A,b);
	Vector<4> u = LT*b;

    this->initial_state[3] += stepSize*u[0];
    this->initial_state[4] += stepSize*u[1];
    this->initial_state[5] += stepSize*u[2];
    this->t_final += stepSize*0.5*u[3];	 

	return delta;
}

void Section::compute_STM(){		
    Matrix<6,6> A;
    Matrix<6,6> dSTM;
    this->STM.set_identity();
    A[0][3] = 1;
    A[1][4] = 1;
    A[2][5] = 1;
    
    const double deltaX = 1;
    
    const std::size_t n = this->times.size() - 1;
    std::array<double,6> a;
    std::array<double,6> x;
    std::array<double,6> ax;
    std::array<double,6> ay;
    std::array<double,6> az;
    for(std::size_t i = 0; i < n;i++) {
        double dt = this->times[i+1] - this->times[i];
        x = this->states[i];
        a = this->dynamics->get_state_rate(x, this->times[i]);
        x[0] += deltaX;
        ax = this->dynamics->get_state_rate(x, this->times[i]);
        x[0] -= deltaX;
        x[1] += deltaX;
        ay = this->dynamics->get_state_rate(x, this->times[i]);
        x[1] -= deltaX;
        x[2] += deltaX;
        az = this->dynamics->get_state_rate(x, this->times[i]);
        
        for(std::size_t j = 3; j < 6; j++){
            ax[j] = (ax[j] - a[j])/deltaX;
            ay[j] = (ay[j] - a[j])/deltaX;
            az[j] = (az[j] - a[j])/deltaX;
        }
        
        for(std::size_t j = 0; j < 3; j++){
            std::size_t j3 = j + 3;
            A[3][j] = ax[j3];
            A[4][j] = ay[j3];
            A[5][j] = az[j3];
        }
        
        dSTM = A*STM;
        STM += (A*STM)*dt;
    }
   
}

void OrbitComputation::differential_correct_cr3bp( CR3BP* cr3bp, std::array<double,6>& x0, std::array<double,6>& xf, double& tf, bool top, const double& tol = 1e-10) {

    double f = cr3bp->sma*cr3bp->mean_motion;
    double tolerance = tol / f;
	double recording_interval = 1e-4; 
	
	ODE_RK4<6> ode = ODE_RK4<6>();
	ode.set_dynamics(cr3bp);
	ode.recording.set_record_interval(recording_interval);
	ode.set_timestep(1e-6);
	if( top ) {
		ode.stop = [](const std::array<double,6>& x,const double& t){
			return x[1] < 0;
		};
	} else {
		ode.stop = [](const std::array<double,6>& x,const double& t){
			return x[1] > 0;
		};
	}
	
	std::cout << "Running" <<std::endl;
	double A00,A11,A10,A01,dx0,dz0,dyd0;
	double** A = Math::zeros(6);
	double** STM = Math::zeros(6);
	double** dSTM = Math::zeros(6);
	bool alternate = true;
	for(int iter = 0; iter < MAX_HALO_ITERATIONS;iter++) {
		ode.recording.clear();
		ode.set_time(0);
		ode.run(x0,100);
		
		xf = ode.get_state();
		
		std::array<double,6> dxf = cr3bp->get_state_rate(xf,ode.get_time());
		
		double dtf = -xf[1]/dxf[1];
		
		double x_dot_f = xf[3] + dxf[3]*dtf;
		double z_dot_f = xf[5] + dxf[5]*dtf;
		
		std::cout << "half period x velocity (m/s): " << x_dot_f*f*1e3 << std::endl;
		std::cout << "half period z velocity (m/s): " << z_dot_f*f*1e3 << std::endl;
		
		double e = (fabs(x_dot_f) + fabs(z_dot_f));
        
        if (e < tolerance) {
			tf = ode.get_time() + dtf;
            break;
		}
    
        //Integrate State Transition Matrix
        Math::identity(STM,6);
        int n = ode.recording.number_entries()-1;
        for(int i = 0; i < n;i++) {
			const std::array<double,6>& xi = ode.recording.state_at(i);
			double dt = ode.recording.time_at(i+1) - ode.recording.time_at(i);
            cr3bp->getA(xi,0,A);
            Math::mult(A,STM,dSTM,6);
			Math::mult(dSTM,dt,6);
            Math::add(STM,dSTM,6);
		}	
		cr3bp->getA(ode.recording.state_at(n),0,A);
		Math::mult(A,STM,dSTM,6);
		Math::mult(dSTM,dtf,6);
		Math::add(STM,dSTM,6);	
    
        //Solve for delta intial state
		
        double xdd = dxf[3]/xf[4];
        double zdd = dxf[5]/xf[4];
    
        //alternate between dx0 and dx0 matrices
        if (alternate) {
            A00 = STM[3][0] - xdd*STM[1][0];
            A01 = STM[3][4] - xdd*STM[1][4];
            A10 = STM[5][0] - zdd*STM[1][0];
            A11 = STM[5][4] - zdd*STM[1][4];
            double det = A00*A11 - A10*A01;
            dx0 = (A01*z_dot_f - A11*x_dot_f)/det;
            dyd0 = (A10*x_dot_f - A00*z_dot_f)/det;
            x0[0] += dx0;
        } else{
            A00 = STM[3][2] - xdd*STM[1][2];
            A01 = STM[3][4] - xdd*STM[1][4];
            A10 = STM[5][2] - zdd*STM[1][2];
            A11 = STM[5][4] - zdd*STM[1][4];
            double det = A00*A11 - A10*A01;
            dz0 = (A01*z_dot_f - A11*x_dot_f)/det;
            dyd0 = (A10*x_dot_f - A00*z_dot_f)/det;
            x0[2] += dz0;
		}
        x0[4] += dyd0;
    
        std::cout << "corrected initial state: " << std::endl;
		
		for(int i = 0; i < 6;i++){
			std::cout << x0[i] << " ";
		}
		std::cout << std::endl;
    
        alternate = !alternate;
	}
	
	Math::del(A,6);
	Math::del(STM,6);
	Math::del(dSTM,6);
}

Recording<6> OrbitComputation::get_cr3bp_halo_orbit(CR3BP* cr3bp, std::array<double,6> x, double dt, double t_period){
	dt *= cr3bp->mean_motion; 
    t_period *= 0.9*cr3bp->mean_motion;
	
	ODE_RK4<6> ode = ODE_RK4<6>();
	ode.set_dynamics(cr3bp);
	ode.recording.set_record_interval(dt);
	ode.set_timestep(1e-6);
	ode.stop = [t_period](const std::array<double,6>& x,const double& t){
		return x[1] > 0 && t > t_period;
	};  

    ode.recording.clear();
    ode.run(x,100);
    
    return ode.recording;
}


void OrbitComputation::minimizeDX(std::vector<Section>& sections){
    std::cout << "patching...\n";
	const double threshold = 1;
    const uint_fast16_t nSections = sections.size(); 
	std::vector<bool> completed(nSections,false);
	int completedCount = 0;
    for(int iter = 0; iter < 30; iter++){
        std::cout << "Iteration: " << iter << std::endl;

        for (uint_fast16_t sId = 0; sId < nSections - 1; sId++) {

			if(completed[sId]){
				continue;
			}

            double dist = sections[sId].target(sections[sId+1].initial_state);

			if(dist < threshold){
				completed[sId] = true;
				completedCount++;
			}

			sections[sId+1].t_start = sections[sId].t_final;
        }

		std::cout << "Completed: " << completedCount << std::endl;

		if(completedCount == nSections-1) {
			break;
		}
    }

	std::cout << "patched.\n";
}

void OrbitComputation::minimizeDV(std::vector<Section>& sections){
    const uint_fast16_t nSections = sections.size();
	std::cout << "Minimizing..." << std::endl;
	
	// Recompute STMS, might not be necessary
	EarthMoonSun* dynamics = sections[0].dynamics;
	for (uint_fast16_t section = 0; section < nSections; section++) {
		sections[section].compute_states();
		sections[section].compute_STM();
	}
	
	Matrix<3,3> Apf;
	Matrix<3,3> Bpf;
	Matrix<3,3> Ap0;
	Matrix<3,3> Bp0;
	Matrix<3,3> M0;
	Vector<3> Mt0;
	Matrix<3,3> Mp;
	Vector<3> Mtp;
	Matrix<3,3> Mf;
	Vector<3> Mtf;
	
	Vector<3> v_back;
	Vector<3> vp_back;
	Vector<3> vp_front;
	Vector<3> ap_back;
	Vector<3> ap_front;
	Vector<3> v_front;
	
	MatrixX M(3*nSections-6, 4*nSections + 4,(char)0);
	MatrixX dvp(M.nRows,1);
	MatrixX dr(M.nCols,1);
		
	for (uint_fast16_t section = 1; section < nSections-1; section++) {
		uint_fast16_t idx_back = section-1;

		uint_fast16_t colStart = idx_back*4;
		uint_fast16_t rowStart = idx_back*3;

		Matrix<6,6> STM_back(sections[idx_back].STM);
		
		Ap0 = STM_back.slice<3,3>(0,0);
		Bp0 = MatrixX::invert(STM_back.slice<3,3>(0,3));
		
		Matrix<6,6> STM_front(sections[section].STM);
		STM_front = MatrixX::LUPInvert(STM_front);
		
		Apf = STM_front.slice<3,3>(0,0);
		Bpf = MatrixX::invert(STM_front.slice<3,3>(0,3));

		const std::array<double,6>& x_back = sections[idx_back].initial_state;
		const std::array<double,6>& xp_back  = dynamics->get_state_rate(sections[idx_back].states.back(), sections[idx_back].times.back());
		const std::array<double,6>& xp_front = dynamics->get_state_rate(sections[section].initial_state, sections[section].t_start);
		const std::array<double,6>& x_front = sections[section].states.back();

		dvp.data[rowStart] = (xp_front[0] - xp_back[0]);
		dvp.data[rowStart+1] = (xp_front[1] - xp_back[1]);
		dvp.data[rowStart+2] = (xp_front[2] - xp_back[2]);

		for(uint_fast8_t i = 0; i < 3; i++){
			v_back[i] = x_back[i+3];
			vp_back[i] = xp_back[i];
			vp_front[i] = xp_front[i];
			v_front[i] = x_front[i+3];
			ap_back[i] = xp_back[i+3];
			ap_front[i] = xp_front[i+3];
		}

		const Matrix<3,3>& tmp0 = Bp0*Ap0;
		const Matrix<3,3>& tmpf = Bpf*Apf;
		
		Mt0 = Bp0*v_back;

		Mtf = Bpf*v_front;

		for(uint_fast8_t i = 0; i < 9; i++){
			Mp.data[i] = tmp0.data[i] - tmpf.data[i];
		}

		const Vector<3>& v0 = tmp0*vp_back;
		const Vector<3>& vf = tmpf*vp_front;

		for(uint_fast8_t i = 0; i < 3; i++){
			Mtp.data[i] = vf.data[i] - v0.data[i] + ap_front.data[i] - ap_back.data[i];				
		}

		for(uint_fast8_t i = 0; i < 3; i++){
			double* row = M[rowStart + i];
			for(uint_fast8_t j = 0; j < 3; j++){
				row[colStart + j] = -Bp0[i][j];
				row[colStart + 4 + j] = Mp[i][j];
				row[colStart + 8 + j] = Bpf[i][j];
			}
			row[colStart + 3] = Mt0(i);
			row[colStart + 7] = Mtp(i);
			row[colStart + 11] = -Mtf(i);
		}
	}
	
	//std::cout << "DV: " << std::endl << dv.to_string() << std::endl;
	
	MatrixX MT = M.get_transpose();
	MatrixX MMT = M*MT;
	
	Math::LUPSolve(MMT.rows,dvp.data,MMT.nRows);

	for(std::size_t i = 0; i < dr.nRows; i++){
		dr.data[i] = Math::dot(MT[i],dvp.data, MT.nCols);
	}

	for (uint_fast16_t section = 0; section < nSections; section++) {
		uint_fast16_t col = section*4;
		for(uint_fast8_t i = 0; i < 3; i++){
			sections[section].initial_state[i] -= 0.5*dr.data[col + i];
		}
		sections[section].t_start -= 0.1*dr.data[col + 3];
	}
	// double check times
	for (uint_fast16_t section = 0; section < nSections-1; section++) {
		sections[section].t_final = sections[section+1].t_start;
	}	

	for (uint_fast16_t section = 0; section < nSections; section++) {
		sections[section].compute_states();
		sections[section].compute_STM();
	}
}


void OrbitComputation::minimizeDV2(std::vector<Section>& sections){
    const uint_fast16_t nSections = sections.size();
	std::cout << "Minimizing..." << std::endl;
	
	

	Matrix<3,3> dVdR0;
	Vector<3> dVdT0;
	Matrix<3,3> dVdV0;
	Vector<3> dVdT;
	
	Vector<3> vp_back;
	Vector<3> vp_front;
	
	Eigen::MatrixXd M(3*nSections-3, 7*nSections);
	Eigen::VectorXd dvp(3*nSections-3);
	Eigen::VectorXd dr(7*nSections);

	M.setZero();

	const double DX = 1;
	const double DV = 0.0001;
	const double DT = 60;
		
	for (uint_fast16_t section = 1; section < nSections; section++) {
		uint_fast16_t idx_back = section-1;

		uint_fast16_t colStart = idx_back*7;
		uint_fast16_t rowStart = idx_back*3;

		std::array<double,6> xp_back  = sections[idx_back].states.back();
		std::array<double,6> xp_front = sections[section].initial_state;

		dvp(rowStart) = (xp_front[3] - xp_back[3]);
		dvp(rowStart+1) = (xp_front[4] - xp_back[4]);
		dvp(rowStart+2) = (xp_front[5] - xp_back[6]);

		for(int i = 0; i < 3; i++){
			sections[idx_back].initial_state[i] += DX;
			sections[idx_back].compute_states();
			sections[idx_back].initial_state[i] -= DX;
			const std::array<double,6>& x_back  = sections[idx_back].states.back();
			for(int j = 0; j < 3;j++){
				dVdR0[j][i] = (x_back[3+j] - xp_back[3+j])/DX;
			}
		}

		for(int i = 0; i < 3; i++){
			sections[idx_back].initial_state[3+i] += DV;
			sections[idx_back].compute_states();
			sections[idx_back].initial_state[3+i] -= DV;
			const std::array<double,6>& x_back  = sections[idx_back].states.back();
			for(int j = 0; j < 3;j++){
				dVdV0[j][i] = (x_back[3+j] - xp_back[3+j])/DV;
			}
		}

		sections[idx_back].t_start += DT;
		sections[idx_back].compute_states();
		sections[idx_back].t_start -= DT;
		const std::array<double,6>& x_back  = sections[idx_back].states.back();
		for(int j = 0; j < 3;j++){
			dVdT0.data[j] = (x_back[3+j] - xp_back[3+j])/DT;
		}

		sections[idx_back].t_final += DT;
		sections[idx_back].compute_states();
		sections[idx_back].t_final -= DT;
		const std::array<double,6>& x_back2  = sections[idx_back].states.back();
		for(int j = 0; j < 3;j++){
			dVdT.data[j] = (x_back2[3+j] - xp_back[3+j])/DT;
		}

		for(uint_fast8_t i = 0; i < 3; i++){
			int row = rowStart + i;
			for(uint_fast8_t j = 0; j < 3; j++){
				M(row,colStart + j) = dVdR0[i][j];
				M(row,colStart + 3 + j) = dVdV0[i][j];
			}
			M(row,colStart + 10 + i) = 1;
			M(row,colStart + 6) = dVdT0(i);
			M(row,colStart + 13) = dVdT(i);
		}
	}

	dr = M.colPivHouseholderQr().solve(dvp);

	for (uint_fast16_t section = 0; section < nSections; section++) {
		uint_fast16_t col = section*7;
		for(uint_fast8_t i = 0; i < 6; i++){
			sections[section].initial_state[i] += 1*dr(col + i);
		}
		sections[section].t_start += 0.5*dr(col + 6);
	}

	for(int i = 0; i < nSections*7;i++){
		std::cout << dr(i) << std::endl;
	}

	// double check times
	for (uint_fast16_t section = 0; section < nSections-1; section++) {
		sections[section].t_final = sections[section+1].t_start;
	}

	for (uint_fast16_t section = 0; section < nSections; section++) {
		sections[section].compute_states();
	}	
}

void OrbitComputation::smooth(std::vector<Section>& sections){
    
}

double OrbitComputation::calcDV(std::vector<Section>& sections){
    double totalDV = 0;
    const uint_fast16_t nSections = sections.size() - 1; 
	for (uint_fast16_t section = 0; section < nSections; section++) {
		std::array<double,6> x_prev = sections[section].states.back();
        std::array<double,6> x_next = sections[section+1].initial_state;
        double dv2 = 0;
        for(uint_fast8_t i = 3; i < 6; i++){
            double dx = x_prev[i] - x_next[i];
            dv2 += dx*dx;
        }
        totalDV += sqrt(dv2);
	}
	return totalDV;
}


void OrbitComputation::run_full_emphemeris(const int& nSections, const double& jd0) {
	EarthMoonSun* dynamics = new EarthMoonSun(jd0);
	
	std::cout << std::setprecision(12);

	// init sections
	std::cout << "Initializing Segments" << std::endl;

    //CR3BP cr3bp = CR3BP(OrbitalElements::EARTH_MU,OrbitalElements::MOON_MU,sma);

	double T = 0;
	const double dT = Section::SECTION_DAYS*86400;
	std::vector<Section> sections(nSections,dynamics);


	for (int section = 0; section < nSections; section++) { 
	
		double jd = jd0 + T/86400.0;
	
		std::array<double,3> e = dynamics->earth->getPos(jd);
		std::array<double,3> m = dynamics->moon->getPos(jd);
		std::array<double,3> r = {e[0] - m[0],e[1] - m[1],e[2] - m[2]};
		double sma = sqrt(Math::dot(r,r));

		
		CR3BP cr3bp = CR3BP(OrbitalElements::EARTH_MU,OrbitalElements::MOON_MU,sma);
		sections[section].initial_state = OrbitalDynamics::convert_cr3bp_to_inertial(dynamics, cr3bp.get_halo_initial_state_3rd_order(10000,0,T,1),jd);
		sections[section].t_start = T;
		T += dT;
		sections[section].t_final = T;
	}
	
	double t_final = T + 3600;
	
	double time = 0;
	std::vector<double> t;
	std::vector<std::array<double,3> > xE;
	std::vector<std::array<double,3> > xM;
	while(time < t_final){
		double jd_t = jd0 + time/86400;
		t.push_back(time);
		xE.push_back(dynamics->earth->getPos(jd_t));
		xM.push_back(dynamics->moon->getPos(jd_t));
		time += 3600;
	}

	std::cout << "Printing OG" << std::endl;
	
	Util::printOut(t,xE,"output/test_earth.orbit");
	Util::printOut(t,xM,"output/test_moon.orbit");

	for (Section& section : sections) { 
		section.compute_states();
	}
	//Util::printOut(dynamics,sections,"output/test_orbit_init.orbit");

	double dvEnd = nSections*0.01;
	const int MAX_ITER = 0;
	for(int iter = 0; iter < MAX_ITER; iter++) {

		OrbitComputation::minimizeDX(sections);
		
		OrbitComputation::minimizeDV(sections);
		
		double dv = OrbitComputation::calcDV(sections);
		std::cout << "DV: " << dv << std::endl;
		if(dv < dvEnd){
			break;
		}

	}

	OrbitComputation::minimizeDX(sections);
	
	std::cout << "printing" << std::endl;
	//Util::printOut(dynamics,sections,"output/test_orbit.orbit");

}

void OrbitComputation::correct(std::vector<Section>& sections){
	for (Section& section : sections) { 
		section.compute_states();
	}
	//Util::printOut(dynamics,sections,"output/test_orbit_init.orbit");

	double dvEnd = sections.size()*0.01;
	const int MAX_ITER = 0;
	for(int iter = 0; iter < MAX_ITER; iter++) {

		OrbitComputation::minimizeDX(sections);
		
		OrbitComputation::minimizeDV(sections);
		
		double dv = OrbitComputation::calcDV(sections);
		std::cout << "DV: " << dv << std::endl;
		if(dv < dvEnd){
			break;
		}

	}

	OrbitComputation::minimizeDX(sections);
}



/*
void compute_STM(){
    Matrix<6,6> A;
    this->STM.set_identity();
    const int n = this->times.size()-1;
    for(int i = 0; i < n;i++) {
        double dt = this->times[i+1] - this->times[i];
        double jd = this->dynamics->getJD0() + this->times[i]/86400;
        this->dynamics->getA(this->states[i],jd,A);
        STM += A*STM*dt;
    }
}
 */
