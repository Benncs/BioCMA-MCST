// #ifndef __MC_PARTICLES_LIST_HPP__
// #define __MC_PARTICLES_LIST_HPP__

// #include <span>
// #include <vector>

// #include <mc/particles/mcparticles.hpp>

// #define __MC_MAX_PARTICLES_PER_PROCESS__ 100e6

// namespace MC
// {

//   static constexpr size_t MAX_PARTICLE_BUFFER =
//       __MC_MAX_PARTICLES_PER_PROCESS__;

//   class ParticlesList
//   {
//   public:
//     explicit ParticlesList() = default;

//     explicit ParticlesList(size_t capacity, double weight);

//     explicit ParticlesList(const ParticlesList &other) = delete;
//     ParticlesList(ParticlesList &&other) noexcept;
//     ~ParticlesList() = default;

//     /*std::vector forward */
//     void emplace_back(Particle &&p);

//     void insert(std::vector<MC::Particle> &&source) noexcept;

//     // auto data(){return data;}

//     template <typename IT> void insert(IT &&begin, IT &&end) noexcept
//     {
//       m_data.insert(m_data.end(),
//                   std::make_move_iterator(begin),
//                   std::make_move_iterator(end));
//     }

//     auto data(){return m_data;}

//     std::span<Particles> data_span(){return m_data;}

//     [[nodiscard]] size_t size() const noexcept;

//     [[nodiscard]] decltype(auto) begin() const
//     {
//       return m_data.begin();
//     }
//     [[nodiscard]] decltype(auto) end() const
//     {
//       return m_data.end();
//     }

//     decltype(auto) begin()
//     {
//       return m_data.begin();
//     }
//     decltype(auto) end()
//     {
//       return m_data.end();
//     }

//     auto &operator[](size_t i);

//     auto &operator[](size_t i) const;

//     ParticlesList &operator=(const ParticlesList &other) = delete;

//     ParticlesList &operator=(ParticlesList &&other) noexcept;

//     template <class Archive> void serialize(Archive &ar) { ar(m_data); }

//   private:
//     std::vector<Particles> m_data;
//   };

//   inline void ParticlesList::emplace_back(Particles &&p)
//   {
//     this->m_data.emplace_back(std::move(p));
//   }

//   inline size_t ParticlesList::size() const noexcept
//   {
//     return m_data.size();
//   }

//   inline auto &ParticlesList::operator[](size_t i)
//   {
//     return m_data[i];
//   }

//   inline auto &ParticlesList::operator[](size_t i) const
//   {
//     return m_data[i];
//   }

//   inline void
//   ParticlesList::insert(std::vector<MC::Particles> &&source) noexcept
//   {
//     m_data.insert(m_data.end(),
//                 std::make_move_iterator(source.begin()),
//                 std::make_move_iterator(source.end()));
//   }
// } // namespace MC

// #endif