#ifndef MYTENSOR_HPP_
#define MYTENSOR_HPP_

#include <iostream>
#include <vector>
#include <functional>
#include <numeric>
#include <cmath>

static bool default_check = false; // set true when you are developing a model

template <class T>
class MyTensor{
private:
    std::vector<size_t> dim;
	std::vector<T> data;
    bool check;
    void error_message();
    T get_data_ind(const std::vector<size_t> arg_ind);

public:
	MyTensor(); // default constructor
	MyTensor(const std::vector<size_t>& arg_dim, const T& c=0.0, const bool& arg_check=default_check);
	MyTensor(size_t i, const T& c=0.0, bool arg_check=default_check);
	MyTensor(size_t i, size_t j, const T& c=0.0, bool arg_check=default_check);
	MyTensor(size_t i, size_t j, size_t k, const T& c=0.0, bool arg_check=default_check);
	MyTensor(size_t i, size_t j, size_t k, size_t l, const T& c=0.0, bool arg_check=default_check);
	MyTensor(size_t i, size_t j, size_t k, size_t l, size_t m, const T& c=0.0, bool arg_check=default_check);
	void resize(const std::vector<size_t>& arg_dim, const T& c=0.0, const bool& arg_check=default_check);
    T& get(const std::vector<size_t>& arg_dim);
    T& get(size_t i);
    T& get(size_t i,size_t j);
    T& get(size_t i,size_t j,size_t k);
    T& get(size_t i,size_t j,size_t k,size_t l);
    T& get(size_t i,size_t j,size_t k,size_t l,size_t m);
    T mean(size_t i, size_t j, size_t k, size_t ii, size_t jj, size_t kk);
    T std(size_t i, size_t j, size_t k, size_t ii, size_t jj, size_t kk);
    T mean(const std::vector<size_t>& sv, const std::vector<size_t>& ev);
    T std(const std::vector<size_t>& sv, const std::vector<size_t>& ev);
    std::vector<size_t>& size();
    MyTensor<T> get_cholesky_lower();
    MyTensor<T>& operator=(const MyTensor<T>& other){
    if (this != &other){ // not a self-assignment
        this->dim = other.dim;
        this->data = other.data;
        this->check = other.check;
        }
        return *this;
    }
    MyTensor<T>& operator=(const T& other){
        for (size_t i=0; i<this->data.size(); i++){
            data[i] = other;
        }
        return *this;
    }
    void print();
    virtual ~MyTensor(){};
};

// Private methods
template <class T>
T MyTensor<T>::get_data_ind(const std::vector<size_t> arg_ind){
    size_t data_ind = 0;
    size_t multiplier = 1;
    for (size_t pntr = arg_ind.size(); pntr --> 0;){
        if (check){
            if (arg_ind[pntr]>=this->dim[pntr]){
                this->error_message();
            }
        }
        data_ind += multiplier*arg_ind[pntr];
        multiplier *= this->dim[pntr];
    }    
    return data_ind;
}

// Public methods

template <class T>
MyTensor<T>::MyTensor(){
    dim.resize(0);
    check = true;
    data.resize(0);
}

template <class T>
MyTensor<T>::MyTensor(const std::vector<size_t>& arg_dim, const T& c, const bool& arg_check){
    dim = arg_dim;
    check = arg_check;
    size_t multi = std::accumulate(dim.begin(), dim.end(), 1, std::multiplies<size_t>());
    data.resize(multi, c);
}

template <class T>
void MyTensor<T>::resize(const std::vector<size_t>& arg_dim, const T& c, const bool& arg_check){
    dim = arg_dim;
    check = arg_check;
    size_t multi = std::accumulate(dim.begin(), dim.end(), 1, std::multiplies<size_t>());
    data.resize(multi, c);
}

template <class T>
MyTensor<T>::MyTensor(size_t i, const T& c, bool arg_check){
    dim.emplace_back(i);
    check = arg_check;
    size_t multi = std::accumulate(dim.begin(), dim.end(), 1, std::multiplies<size_t>());
    data.resize(multi, c);
}

template <class T>
MyTensor<T>::MyTensor(size_t i, size_t j, const T& c, bool arg_check){
    dim.emplace_back(i);
    dim.emplace_back(j);
    check = arg_check;
    size_t multi = std::accumulate(dim.begin(), dim.end(), 1, std::multiplies<size_t>());
    data.resize(multi, c);
}
	
template <class T>
MyTensor<T>::MyTensor(size_t i, size_t j, size_t k, const T& c, bool arg_check){
    dim.emplace_back(i);
    dim.emplace_back(j);
    dim.emplace_back(k);
    check = arg_check;
    size_t multi = std::accumulate(dim.begin(), dim.end(), 1, std::multiplies<size_t>());
    data.resize(multi, c);
}

template <class T>
MyTensor<T>::MyTensor(size_t i, size_t j, size_t k, size_t l, const T& c, bool arg_check){
    dim.emplace_back(i);
    dim.emplace_back(j);
    dim.emplace_back(k);
    dim.emplace_back(l);
    check = arg_check;
    size_t multi = std::accumulate(dim.begin(), dim.end(), 1, std::multiplies<size_t>());
    data.resize(multi, c);
}

template <class T>
MyTensor<T>::MyTensor(size_t i, size_t j, size_t k, size_t l, size_t m, const T& c, bool arg_check){
    dim.emplace_back(i);
    dim.emplace_back(j);
    dim.emplace_back(k);
    dim.emplace_back(l);
    dim.emplace_back(m);
    check = arg_check;
    size_t multi = std::accumulate(dim.begin(), dim.end(), 1, std::multiplies<size_t>());
    data.resize(multi, c);
}

template <class T>
void MyTensor<T>::error_message(){
    std::cout << "Violation of dimension in MyTensor!" << std::endl;
    exit(1);
}

template <class T>
T& MyTensor<T>::get(const std::vector<size_t>& arg_ind){
    return data[get_data_ind(arg_ind)];
}
    
template <class T> 
T& MyTensor<T>::get(size_t i){
    if (check){
        if (i>=dim[0]){
            error_message();
        }
    }
    return data[i];
}

template <class T> 
T& MyTensor<T>::get(size_t i,size_t j){
    if (check){
        if ((i>=dim[0])||(j>=dim[1])){
            error_message();
        }
    }
    return data[i*dim[1]+j];
}

template <class T> 
T& MyTensor<T>::get(size_t i,size_t j,size_t k){
    if (check){
        if ((i>=dim[0])||(j>=dim[1])||(k>=dim[2])){
            error_message();
        }
    }
    return data[i*dim[1]*dim[2]+j*dim[2]+k];
}

template <class T> 
T& MyTensor<T>::get(size_t i,size_t j,size_t k,size_t l){
    if (check){
        if ((i>=dim[0])||(j>=dim[1])||(k>=dim[2])||(l=dim[3])){
            error_message();
        }
    }    
    return data[i*dim[1]*dim[2]*dim[3]+j*dim[2]*dim[3]+k*dim[3]+l];
}

template <class T> 
T& MyTensor<T>::get(size_t i,size_t j,size_t k,size_t l,size_t m){
    if (check){
        if ((i>=dim[0])||(j>=dim[1])||(k>=dim[2])||(l=dim[3])||(m=dim[4])){
            error_message();
        }
    }
    return data[i*dim[1]*dim[2]*dim[3]*dim[4]+j*dim[2]*dim[3]*dim[4]+k*dim[3]*dim[4]+l*dim[4]+m];
}

template <class T> 
MyTensor<T> MyTensor<T>::get_cholesky_lower(){
    MyTensor lower = MyTensor(dim, static_cast<T>(0.0)); // Initialize the result lower matrix
    for (size_t i=0; i<dim[0]; i++){
        for (size_t j=0; j<=i; j++){
            T sum = 0;
            if (j==i){ // summation for diagonals
                for (size_t k = 0; k < j; k++){
                    sum += std::pow(lower.get(j,k), 2);
                } 
                lower.get(j,j) = std::sqrt(this->get(j,j) - sum);
            }
            else {
                // Evaluating L(i, j) using L(j, j)
                for (size_t k = 0; k < j; k++) {
                    sum += (lower.get(i,k) * lower.get(j,k));
                }
                lower.get(i,j) = (this->get(i,j) - sum) / lower.get(j,j);
            }
        }
    }
    return lower;

    // Cholesky test
    // MyTensor test = MyTensor(static_cast<size_t>(3),static_cast<size_t>(3),0.0);
    // test.get(0,0)=4.0; test.get(0,1)=12.0; test.get(0,2)=-16.0;
    // test.get(1,0)=12.0; test.get(1,1)=37.0; test.get(1,2)=-43.0;
    // test.get(2,0)=-16.0; test.get(2,1)=-43.0; test.get(2,2)=98.0;
    // test.print(); std::cout << std::endl;
}

template <class T> 
void MyTensor<T>::print(){
    if (dim.size()==1){
        for(size_t i=0; i<dim[0]; i++){
            std::cout << data[i] << std::endl;
        }
    }
    else if (dim.size()==2){
        for(size_t i=0; i<dim[0]; i++){
            for(size_t j=0; j<dim[1]; j++){
                std::cout << data[i*dim[1]+j] << " ";
            }
            std::cout << std::endl;
        }
    }
}

template <class T> 
std::vector<size_t>& MyTensor<T>::size(){
    return dim;
}

template <class T> 
T MyTensor<T>::mean(size_t i, size_t j, size_t k, size_t ii, size_t jj, size_t kk){
    size_t s = i*dim[1]*dim[2] + j*dim[2] + k;
    size_t e = ii*dim[1]*dim[2] + jj*dim[2] + kk;
    if (e <= s){
        std::cout << "Wrong index for MyTensor<T>::mean" << std::endl;
        exit(-1);
    }
    T sum = std::accumulate(&data[s], &data[e+1], static_cast<T>(0.0));
    T mean = sum / static_cast<T>(e-s+1);
    // std::cout << data[s] << " " << data[e]  << std::endl;
    // std::cout << sum << " " << mean << " " << s << " " << e << std::endl;
    return mean;
}

template <class T> 
T MyTensor<T>::std(size_t i, size_t j, size_t k, size_t ii, size_t jj, size_t kk){
    size_t s = i*dim[1]*dim[2] + j*dim[2] + k;
    size_t e = ii*dim[1]*dim[2] + jj*dim[2] + kk;
    if (e <= s){
        std::cout << "Wrong index for MyTensor<T>::mean" << std::endl;
        exit(-1);
    }
    T m = this->mean(i,j,k,ii,jj,kk);

    T sq_sum = std::inner_product(&data[s], &data[e+1], &data[s], static_cast<T>(0.0));
    T stdev = std::sqrt(sq_sum/static_cast<T>(e-s+1) - m*m);

    // std::cout << data[s] << " " << data[e]  << std::endl;
    // std::cout << sum << " " << mean << " " << s << " " << e << std::endl;
    return stdev;
}

template <class T>
T MyTensor<T>::mean(const std::vector<size_t>& sv, const std::vector<size_t>& ev){
    size_t s = get_data_ind(sv);
    size_t e = get_data_ind(ev);
    if (e <= s){
        std::cout << "Wrong index for MyTensor<T>::mean" << std::endl;
        exit(-1);
    }
    T sum = std::accumulate(&data[s], &data[e+1], static_cast<T>(0.0));
    T mean = sum / static_cast<T>(e-s+1);
    // // std::cout << data[s] << " " << data[e]  << std::endl;
    // // std::cout << sum << " " << mean << " " << s << " " << e << std::endl;
    return mean;
}

template <class T>
T MyTensor<T>::std(const std::vector<size_t>& sv, const std::vector<size_t>& ev){
    size_t s = get_data_ind(sv);
    size_t e = get_data_ind(ev);
    if (e <= s){
        std::cout << "Wrong index for MyTensor<T>::mean" << std::endl;
        exit(-1);
    }
    T m = this->mean(sv, ev);

    T sq_sum = std::inner_product(&data[s], &data[e+1], &data[s], static_cast<T>(0.0));
    T stdev = std::sqrt(sq_sum/static_cast<T>(e-s+1) - m*m);

    // std::cout << data[s] << " " << data[e]  << std::endl;
    // std::cout << sum << " " << mean << " " << s << " " << e << std::endl;
    return stdev;
}


#endif /* MYTENSOR_HPP_ */