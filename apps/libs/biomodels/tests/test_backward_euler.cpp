
#include <backward_euler.hpp>
#include <Eigen/Core>
#include <iostream>
#include <cassert> 


namespace {
    using xi_t = Eigen::Vector<double, 2>;
}

struct DummyModel
{
    xi_t xi;
};


void f(double s, const xi_t &xi, xi_t &xi_dot)
{

    xi_dot = -15*xi;

}

int main()
{
    double t =0;
    double d_t = 0.001;

    DummyModel model;
    model.xi.setOnes();
    size_t n_t = 40; 

    for(size_t i =0;i<n_t;++i)
    {
        xi_t xi_prev = model.xi; 
        model.xi = backward_euler_update<2>(d_t, t, model.xi, f);
        t += d_t;

        // assert((model.xi.array() >= xi_prev.array()).all());
    }

    double expected_value = std::exp(-15*t);
    std::cout<<model.xi[0] <<"\t"<< expected_value<<std::endl;
    assert(std::abs(model.xi[0] - expected_value)/std::abs(expected_value)*100. < 1); //1% error 
    assert(std::abs(model.xi[1] - expected_value)/std::abs(expected_value)*100. < 1); 
}