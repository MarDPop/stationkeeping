#include "../include/Simulation.h"
#include "../include/ODE.h"
#include "../include/Util.h"
#include "../include/OrbitalDynamics.h"
#include "../include/OrbitalElements.h"

#define MAX_HALO_ITERATIONS 10


const double Section::SECTION_DAYS = 1;

Section::Section(EarthMoonSun* f){
    this->dynamics = f;
    this->ode = ODE_RK4<6>();
    this->ode.set_dynamics(f);
    this->ode.recording.set_record_interval(100);
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

void Section::target(const std::array<double,6>& xp){
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

    Matrix<4,3> LT = L.get_transpose();
    Matrix<3,3> A = L*LT;
    Matrix<3,3>::solve(A,b);
	Vector<4> u = LT*b;

    this->initial_state[3] += 0.2*u[0];
    this->initial_state[4] += 0.2*u[1];
    this->initial_state[5] += 0.2*u[2];
    this->t_final += 0.1*u[3];	
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
			const std::array<double,6>& xi = ode.recording.get(i);
			double dt = ode.recording.time_at(i+1) - ode.recording.time_at(i);
            cr3bp->getA(xi,0,A);
            Math::mult(A,STM,dSTM,6);
			Math::mult(dSTM,dt,6);
            Math::add(STM,dSTM,6);
		}	
		cr3bp->getA(ode.recording.get(n),0,A);
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
    const uint_fast16_t nSections = sections.size(); 
    for(int iter = 0; iter < 20; iter++){
        std::cout << "Iteration: " << iter << std::endl;

        for (uint_fast16_t sId = 0; sId < nSections - 1; sId++) {
            sections[sId].target(sections[sId+1].initial_state);
        }

        for (uint_fast16_t sId = 1; sId < nSections; sId++) {
            sections[sId].t_start = sections[sId-1].t_final;
        }
    }
	
	
	std::cout << "patched.\n";
}

void OrbitComputation::minimizeDV(std::vector<Section>& sections){
    const uint_fast16_t nSections = sections.size();
	std::cout << "Minimizing..." << std::endl;
	
	// Recompute STMS, might not be necessary
	/*
	for (uint_fast16_t section = 0; section < nSections; section++) {
		sections[section].compute_states();
		sections[section].compute_STM();
	}
	*/
	
	Matrix<3,3> Apf;
	Matrix<3,3> Bpf;
	Matrix<3,3> Cpf;
	Matrix<3,3> Dpf;
	Matrix<3,3> Ap0;
	Matrix<3,3> Bp0;
	Matrix<3,3> Cp0;
	Matrix<3,3> Dp0;
	Matrix<3,3> M0;
	Vector<3> Mt0;
	Matrix<3,3> Mp;
	Vector<3> Mtp;
	Matrix<3,3> Mf;
	Vector<3> Mtf;
	
	Vector<3> vp_back;
	Vector<3> vp_front;
	Vector<3> ap_back;
	Vector<3> ap_front;
	
	MatrixX M(3*nSections-3, 4*nSections + 4,(char)0);
	MatrixX dvp(M.nRows,1);
	MatrixX dr(M.nCols,1);
	MatrixX dv(nSections-1,1);
		
	for (uint_fast16_t section = 0; section < nSections-1; section++) {
		uint_fast16_t colStart = section*4;
		uint_fast16_t rowStart = section*3;
		
		Ap0 = sections[section].STM.slice<3,3>(0,0);
		Bp0 = sections[section].STM.slice<3,3>(0,3);
		Cp0 = sections[section].STM.slice<3,3>(3,0);
		Dp0 = sections[section].STM.slice<3,3>(3,3);
		
		Matrix<6,6> STM_front(sections[section+1].STM);
		Matrix<6,6> STM_reverse = MatrixX::LUPInvert(STM_front);
		
		Apf = STM_reverse.slice<3,3>(0,0);
		Bpf = STM_reverse.slice<3,3>(0,3);
		Cpf = STM_reverse.slice<3,3>(3,0);
		Dpf = STM_reverse.slice<3,3>(3,3);

		const std::array<double,6>& x_back  = sections[section].dynamics->get_state_rate(sections[section].states.back(), sections[section].times.back());
		const std::array<double,6>& x_front = sections[section+1].dynamics->get_state_rate(sections[section+1].initial_state, sections[section+1].t_start);
		
		dvp.data[rowStart] = (x_front[0] - x_back[0]);
		dvp.data[rowStart+1] = (x_front[1] - x_back[1]);
		dvp.data[rowStart+2] = (x_front[2] - x_back[2]);

		for(uint_fast8_t i = 0; i < 3; i++){
			vp_back[i] = x_back[i];
			vp_front[i] = x_front[i];
			ap_back[i] = x_back[i+3];
			ap_front[i] = x_front[i+3];
		}

		const Matrix<3,3>& tmp0 = Dp0*MatrixX::invert(Bp0);
		const Matrix<3,3>& tmpf = Dpf*MatrixX::invert(Bpf);
		
		M0 = tmp0*Ap0 - Cp0;
		Mt0 = ap_back - tmp0*vp_back;

		Mf = Cpf - tmpf*Apf;
		Mtf = tmpf*vp_front - ap_front;
		
		for(uint_fast8_t i = 0; i < 3; i++){
			Mtp.data[i] = - Mt0.data[i] - Mtf.data[i];				
		}
		for(uint_fast8_t i = 0; i < 9; i++){
			Mp.data[i] = tmpf.data[i] - tmp0.data[i];
		}

		for(uint_fast8_t i = 0; i < 3; i++){
			double* row = M[rowStart + i];
			for(uint_fast8_t j = 0; j < 3; j++){
				row[colStart + j] = M0[i][j];
				row[colStart + 4 + j] = Mp[i][j];
				row[colStart + 8 + j] = Mf[i][j];
			}
			row[colStart + 3] = Mt0(i);
			row[colStart + 7] = Mtp(i);
			row[colStart + 11] = Mtf(i);
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
			sections[section].initial_state[i] -= dr.data[col + i];
		}
		sections[section].t_start -= 0.5*dr.data[col + 3];
	}
	// double check times
	for (uint_fast16_t section = 0; section < nSections-1; section++) {
		sections[section].t_final = sections[section+1].t_start;
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
		sections[section].initial_state = OrbitComputation::convert(&cr3bp, dynamics, cr3bp.get_halo_initial_state_3rd_order(10000,0,T,1),jd);
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

std::array<double,6> OrbitComputation::convert(CR3BP* cr3bp, EarthMoonSun* dynamics, const std::array<double,6>& state_guess, const double& jd){
	std::array<double,6> x = cr3bp->convert_state_to_inertial(state_guess);
	
	const double xL1 = cr3bp->getL1() - cr3bp->mu;
	
	std::array< std::array<double,3>, 4> frame = dynamics->getEarthMoonBarycenterCS(jd);
	
	std::array<double,3> moon = dynamics->moon->getPos(jd - 0.0005);
	std::array<double,3> earth = dynamics->earth->getPos(jd - 0.0005);
	std::array<double,3> L1_0 = {moon[0] - earth[0],moon[1] - earth[1],moon[2] - earth[2]};
	
	moon = dynamics->moon->getPos(jd + 0.0005);
	earth = dynamics->earth->getPos(jd + 0.0005);
	std::array<double,3> L1_1 = {moon[0] - earth[0],moon[1] - earth[1],moon[2] - earth[2]};
	
	double dt = 0.001*86400;
	double dL1dt = (sqrt(Math::dot(L1_1,L1_1)) - sqrt(Math::dot(L1_0,L1_0)))*xL1/dt;
	
	std::array< std::array<double,3>, 3> CS = {frame[1],frame[2],frame[3]};
	std::array< std::array<double,3>, 3> CST = Math::transpose(CS);
	
	std::array<double,3> pos = {x[0],x[1],x[2]};
	std::array<double,3> vel = {x[3] + dL1dt,x[4],x[5]};
	
	pos = Math::mult(CST,pos);
	vel = Math::mult(CST,vel);

	std::array<double,3>& origin = frame[0]; // Should be zero!
	x[0] = pos[0] + origin[0];
	x[1] = pos[1] + origin[1];
	x[2] = pos[2] + origin[2];
	x[3] = vel[0];
	x[4] = vel[1];
	x[5] = vel[2];
	
	return x;
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
