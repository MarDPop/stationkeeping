#include <vector>
#include <array>

template<uint_fast16_t N>
class DoubleTable {
protected:
	std::vector< double > index;
	std::vector< std::array<double,N> > data;
public:
	
	DoubleTable(){}
	DoubleTable(const std::vector< double >& i, const std::vector< std::array<double,N> >& d) : index(i), data(d) {}
	~DoubleTable(){}
	
	std::vector<double> get_indices() {
		return this->index;
	}
	
	std::vector< std::array<double,N> > get_data() {
		return this->data;
	}
	
	virtual void add(const double& x, const std::array<double,N> value) {
		if(x <= index[0]) {
			index.insert(0,x);
			data.insert(0,value);
			return;
		}
		if(x >= index.back()) {
			index.push_back(x);
			data.push_back(value);
			return;
		}
		uint_fast32_t idx = find(x);
		index.insert(idx,x);
		data.insert(idx,value);
	}
	
	uint_fast32_t find(const double& x) const {
		uint_fast32_t lo = 0;
		uint_fast32_t hi = index.size() - 1;	
		uint_fast32_t mid = (lo + hi)/2;
		while(lo != mid) {
			if (x > index[mid]){
				lo = mid;
			} else {
				hi = mid;
			}
			mid = (lo + hi)/2;
		}
		return mid;
	}
	
	virtual std::array<double,N> get(const double& x) const {
		if(x <= index[0]) {
			return data[0];
		}
		if(x >= index.back()) {
			return data.back();
		}
		
		return data[this->find(x)];
	}
	
}

template<uint_fast16_t N>
class LinearTable : public DoubleTable<N> {
protected:
	std::vector< std::array<double,N> > dvdx;
public:

	

	virtual std::array<double,N> get(const double& x) const {
		if(x <= index[0]) {
			return data[0];
		}
		if(x >= index.back()) {
			return data.back();
		}
		uint_fast32_t = this->find(x);
		std::array<double,N> v;
		for(uint_fast16_t i = 0; i < N; i++){
			
		}
		
		return v;
	}
};

template<uint_fast16_t N>
class CubicTable : public DoubleTable<N> {
protected:
	std::vector< std::array<double,N> > p1;
	std::vector< std::array<double,N> > p2;
	std::vector< std::array<double,N> > p3;
	std::vector< std::array<double,N> > p4;
public:
	
};


template< uint_fast16_t ROW, uint_fast16_t COL >
struct FixedDoubleTable {
	const double index_start;
	const double delta;
	const uint_fast32_t NDATA = ROW*COL;
	double DATA[NDATA] __attribute__ ((aligned (32)));
	
	FixedDoubleTable(const double& d, const double& start) : delta(d), index_start(start) {}
	~FixedDoubleTable(){}
	
	void get(const double& x, double out[COL]) const {
		if(x <= index_start) {
			return data[0];
		}
		if(x >= index.back()) {
			return data.back();
		}
		
		return data[this->find(x)];
	}
	
}