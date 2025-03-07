#include "load_balancing/impl_lb.hpp"
#include <cassert>
#include <iostream>
#include <load_balancing/iload_balancer.hpp>

void test(ILoadBalancer* lb,uint64_t n,uint32_t n_rank)
{
    uint64_t cumsum =0;
    for(uint32_t i=0;i<n_rank;++i)
    {
        auto ib = lb->balance(i, n);
        cumsum+=ib;
        std::cerr<<ib<<std::endl;
    }
    std::cerr<<cumsum<<" "<<n<<std::endl;
    assert(lb->check());
    assert(cumsum==n);
}


int main()
{
    const uint32_t n_rank = 10;
    const uint32_t n_particle = 10e6;
    UniformLoadBalancer uniform(n_rank);
    CustomLoadBalancer custom(n_rank,{0.2,0.15,0.1,0.1,0.2,0.1,0.05,0.05,0.025,0.025});
    HostImportantLoadBalancer host1(n_rank,0.5);
    HostImportantLoadBalancer host2(n_rank,3.1);


    test(&uniform, n_particle , n_rank);
    test(&custom, n_particle , n_rank);
    test(&host1, n_particle , n_rank);
    test(&host2, n_particle , n_rank);

    
}