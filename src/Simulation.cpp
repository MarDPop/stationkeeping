#include "../include/Simulation.h"

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

void Section::minimizeDX(std::vector<Section>& sections){
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

void Section::minimizeDV(std::vector<Section>& sections){
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


void Section::smooth(std::vector<Section>& sections){
    
    for (uint_fast16_t section = 0; section < nSections; section++) {
		sections[section].compute_states();
		sections[section].compute_STM();
	}
}

double Section::calcDV(std::vector<Section>& sections){
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

/*
void compute_STM3(){		
    // Markley Method
    Matrix<6,6> dSTM;
    Matrix<3,3> RR;
    Matrix<3,3> RV;
    Matrix<3,3> VR;
    Matrix<3,3> VV;

    double I[9] = {1,0,0,0,1,0,0,0,1};
    this->STM.set_identity();
    const int n = this->times.size() - 1;
    std::vector< Matrix<3,3> > Gs;
    for (int i = 0; i <= n;i++){
        double jd = this->dynamics->getJD0() + this->times[i]/86400;
        Gs.push_back(std::move(this->dynamics->getG(this->states[i],jd)));
    }	

    for (int i = 0; i < n;i++) {
        double dt = this->times[i+1] - this->times[i];
        double dt_half = dt*0.5;
        double dt6 = dt*dt*0.1666666666666666666666;
        Matrix<3,3>& G0 = Gs[i];
        Matrix<3,3>& G = Gs[i + 1];
        for (int j = 0; j < 9; j++) {
            double tmp = G0.data[j] + G.data[j];
            VR.data[j] = tmp*dt_half;
            RR.data[j] = I[j] + (tmp + G0.data[j])*dt6;
            RV.data[j] = I[j]*dt + VR.data[j]*dt6;
            VV.data[j] = I[j] + (tmp + G.data[j])*dt6; 
        }
        
        for (int j = 3; j < 3; j++) {
            for (int k = 3; k < 3; k++) {
                dSTM.set(i,j,RR(i,j));
                dSTM.set(i,j+3,RV(i,j));
                dSTM.set(i+3,j,VR(i,j));
                dSTM.set(i+3,j+3,VV(i,j));
            }
        }
        
        this->STM = dSTM*this->STM;
    }
}

void compute_STM4(){	
    const double deltaX = 1;
    const double deltaX_half = deltaX*0.5;
    const double deltaV = 0.001;	
    const double deltaV_half = deltaV*0.5;
    
    std::array<double,6> xf0,xf1;
    for(int i = 0; i < 3; i++){
        this->initial_state[i] -= deltaX_half;
        this->compute_states();
        xf0 = this->states.back();
        this->initial_state[i] += deltaX;
        this->compute_states();
        xf1 = this->states.back();
        this->initial_state[i] -= deltaX_half;
        for(int j = 0; j < 6; j++){
            this->STM[j][i] = (xf1[j] - xf0[j])/deltaX;
        }
    }
    
    for(int i = 3; i < 6; i++){
        this->initial_state[i] -= deltaV_half;
        this->compute_states();
        xf0 = this->states.back();
        this->initial_state[i] += deltaV;
        this->compute_states();
        xf1 = this->states.back();
        this->initial_state[i] -= deltaV_half;
        for(int j = 0; j < 6; j++){
            this->STM[j][i] = (xf1[j] - xf0[j])/deltaV;
        }
    }
    this->compute_states();
}
*/