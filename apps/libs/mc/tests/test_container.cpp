#include "Kokkos_Assert.hpp"
#include "mc/traits.hpp"
#include <mc/unit.hpp>
#include <mc/particles_container.hpp>




int main()
{
    const std::size_t size = 1000;
    MC::ParticlesContainer<DefaultModel> container(size);

    KOKKOS_ASSERT(container.n_particles()==size)
    KOKKOS_ASSERT(container.capacity()==size*container.get_allocation_factor());

    
}