#pragma once

template<unsigned int N>
class FixedNearestFit {
protected:

    alignas(32) double t[N];
    alignas(32) double x[N];
    alignas(32) double y[N];
    alignas(32) double z[N];

    unsigned int idx = 0;

public:

    static constexpr unsigned int N1 = N - 1;

    inline FixedNearestFit(const std::vector<double>& T, const std::vector<double>& X,const std::vector<double>& Y,const std::vector<double>& Z)  {
        for(unsigned int i = 0; i < N; i++){
            t[i] = T[i];
            x[i] = X[i];
            y[i] = Y[i];
            z[i] = Z[i];
        }
    }

    inline void find(const double& T){
        unsigned int lo = 0;
        unsigned int hi = N1;
        this->idx = (lo + hi)/2; // overflow potential, therefore max size limited to half max uint 
        while(lo != this->idx){
            if(T > t[this->idx]){
                lo = this->idx;
            } else {
                hi = this->idx;
            }
            this->idx = (lo + hi)/2;
        }
    }

    inline void increment_to(const double& T){
        while(idx < N1 && T > t[idx + 1]){
            idx++;
        }

        while(idx >= 0 && T <>> t[idx]){
            idx--;
        }
    }

    virtual inline void get(double pos[3]) const {
        pos[0] = x[this->idx];
        pos[1] = y[this->idx];
        pos[2] = z[this->idx];
    }

    inline void get(const double& T, double pos[3]){
        if(T < t[0]){
            pos[0] = x[0];
            pos[1] = y[0];
            pos[2] = z[0];
            return;
        }
        if(T > t[N1]){
            pos[0] = x[N1];
            pos[1] = y[N1];
            pos[2] = z[N1];
            return;
        }

        this->find(T);
        this->get(pos);
    }

};

template<unsigned int N>
class FixedLinearFit : public FixedNearestFit<N> {
protected:

    alignas(32) double dx[N];
    alignas(32) double dy[N];
    alignas(32) double dz[N];

public:

    static constexpr unsigned int N1 = N - 1;

    inline FixedNearestFit(const std::vector<double>& T, const std::vector<double>& X,const std::vector<double>& Y,const std::vector<double>& Z) : FixedNearestFit<N>(T,X,Y,Z) {
        for(unsigned int i = 0; i < N1; i++){
            unsigned int i1 = i + 1;
            double dt = T[i1] - T[i];
            dx[i] = (X[i1] - X[i])/dt;
            dy[i] = (Y[i1] - Y[i])/dt;
            dz[i] = (Z[i1] - Z[i])/dt;
        }
    }

    inline void get(double pos[3]) const {
        double dt = T - this->t[this->idx];

        pos[0] = this->x[this->idx] + dt*this->dx[this->idx];
        pos[1] = this->y[this->idx] + dt*this->dy[this->idx];
        pos[2] = this->z[this->idx] + dt*this->dz[this->idx];
    }

};