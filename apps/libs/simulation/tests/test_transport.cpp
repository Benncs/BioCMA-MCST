
#include "simulation/transport.hpp"
#include <Eigen/Dense>
#include <cassert>
#include <vector>

static size_t n_c = 16;
static std::vector<std::vector<double>> get_raw_data();
static std::vector<std::vector<size_t>> get_neighbors();

void check_args();

int main()
{
  check_args();
  auto neighbors = get_neighbors();
  // Initialize an empty vector to store flattened data
  std::vector<double> flattened_data;
  auto data = get_raw_data();
  // Flatten the 2D 'data' vector and store it in 'flattened_data'
  for (const auto &row : data)
  {
    // Append each element of 'row' to 'flattened_data'
    flattened_data.insert(flattened_data.end(), row.begin(), row.end());
  }

  // Map the flattened data to an Eigen Matrix
  Eigen::Map<Eigen::MatrixXd> matrix(
      flattened_data.data(), data.size(), data[0].size());

  // Create a sparse view of the Eigen Matrix
  auto sparse = matrix.sparseView();

  // Create an instance of MatFlow struct
  Simulation::MatFlow f;

  // Set the 'flows' member of 'f' to the sparse view
  f.flows = sparse;

  // Calculate the cumulative probability matrix
  auto cp = Simulation::get_CP(neighbors, n_c, f);

  // Get the maximum index of neighbors
  int neighbor_max = neighbors[0].size() - 1;

  // Loop through each row of the cumulative probability matrix 'cp'
  for (int j = 0; j < n_c; ++j)
  {
    // Initialize the last value to the first element of the row
    double last_value = cp.coeff(j, 0);

    // Assert that the probability is between 0 and 1
    assert(last_value >= 0 && last_value <= 1);

    // Iterate through each neighbor
    for (int k = 1; k < neighbors[0].size(); ++k)
    {
      // Get the current probability value
      double current = cp.coeff(j, k);

      // Assert that the probability increases along the neighbors dimension
      assert(current >= last_value);

      // Assert that the probability is between 0 and 1
      assert(current >= 0 && current <= 1);

      last_value = current;
    }

    // Assert that if there is flow, the last probability is 1; if no flow, it's
    // 0
    assert(cp.coeff(j, neighbor_max) == 1 || cp.coeff(j, neighbor_max) == 0);
  }
}

void check_args()
{
  auto data = get_raw_data();
  auto neighbors = get_neighbors();
  assert(data.size() == n_c);
  for (auto &&i : data)
  {
    assert(i.size() == n_c);
  }

  assert(neighbors.size() == n_c);
}

static std::vector<std::vector<size_t>> get_neighbors()
{
  return {{1, 4, 8},
          {2, 5, 9},
          {3, 6, 10},
          {7, 11, 3},
          {5, 12, 4},
          {6, 13, 5},
          {7, 14, 6},
          {15, 7, 7},
          {9, 12, 8},
          {10, 13, 9},
          {11, 14, 10},
          {15, 11, 11},
          {13, 12, 12},
          {14, 13, 13},
          {15, 14, 14},
          {15, 15, 15}};
}

static std::vector<std::vector<double>> get_raw_data()
{
  return {{0.0,
           0.0,
           0.0,
           0.0,
           0.0132029525,
           0.0,
           0.0,
           0.0,
           3.22623937,
           0.0,
           0.0,
           0.0,
           0.0,
           0.0,
           0.0,
           0.0},
          {0.492171224,
           0.0,
           1.53879051,
           0.0,
           0.0,
           0.00223187291,
           0.0,
           0.0,
           0.0,
           1.04265274,
           0.0,
           0.0,
           0.0,
           0.0,
           0.0,
           0.0},
          {0.0,
           0.171994054,
           0.0,
           1.69650538,
           0.0,
           0.0,
           0.0000025365451,
           0.0,
           0.0,
           0.0,
           1.55723453,
           0.0,
           0.0,
           0.0,
           0.0,
           0.0},
          {0.0,
           0.0,
           0.0492613003,
           0.0,
           0.0,
           0.0,
           0.0,
           0.0136854143,
           0.0,
           0.0,
           0.0,
           2.97005909,
           0.0,
           0.0,
           0.0,
           0.0},
          {4.11125851,
           0.0,
           0.0,
           0.0,
           0.0,
           0.0,
           0.0,
           0.0,
           0.0,
           0.0,
           0.0,
           0.0,
           3.31006412,
           0.0,
           0.0,
           0.0},
          {0.0,
           2.00204127,
           0.0,
           0.0,
           0.773879937,
           0.0,
           1.41600711,
           0.0,
           0.0,
           0.0,
           0.0,
           0.0,
           0.0,
           0.893554508,
           0.0,
           0.0},
          {0.0,
           0.0,
           1.54490607,
           0.0,
           0.0,
           0.119737363,
           0.0,
           1.50000504,
           0.0,
           0.0,
           0.0,
           0.0,
           0.0,
           0.0,
           1.45831214,
           0.0},
          {0.0,
           0.0,
           0.0,
           2.56936589,
           0.0,
           0.0,
           0.0529315415,
           0.0,
           0.0,
           0.0,
           0.0,
           0.0,
           0.0,
           0.0,
           0.0,
           2.6163922},
          {3.28154447,
           0.0,
           0.0,
           0.0,
           0.0,
           0.0,
           0.0,
           0.0,
           0.0,
           0.76409089,
           0.0,
           0.0,
           0.00670166089,
           0.0,
           0.0,
           0.0},
          {0.0,
           2.66697557,
           0.0,
           0.0,
           0.0,
           0.0,
           0.0,
           0.0,
           0.67036601,
           0.0,
           0.418601338,
           0.0,
           0.0,
           0.0519765228,
           0.0,
           0.0},
          {0.0,
           0.0,
           1.82840704,
           0.0,
           0.0,
           0.0,
           0.0,
           0.0,
           0.0,
           1.81558363,
           0.0,
           0.810012534,
           0.0,
           0.0,
           0.0000910195949,
           0.0},
          {0.0,
           0.0,
           0.0,
           1.33402353,
           0.0,
           0.0,
           0.0,
           0.0,
           0.0,
           0.0,
           2.3361757,
           0.0,
           0.0,
           0.0,
           0.0,
           0.00194182384},
          {0.0,
           0.0,
           0.0,
           0.0,
           2.15573001,
           0.0,
           0.0,
           0.0,
           2.10755365,
           0.0,
           0.0,
           0.0,
           0.0,
           1.82208297,
           0.0,
           0.0},
          {0.0,
           0.0,
           0.0,
           0.0,
           0.0,
           2.96898145,
           0.0,
           0.0,
           0.0,
           1.67575778,
           0.0,
           0.0,
           0.780379686,
           0.0,
           0.366231911,
           0.0},
          {0.0,
           0.0,
           0.0,
           0.0,
           0.0,
           0.0,
           1.86931817,
           0.0,
           0.0,
           0.0,
           1.5505214,
           0.0,
           0.0,
           2.00629146,
           0.0,
           0.74915897},
          {0.0,
           0.0,
           0.0,
           0.0,
           0.0,
           0.0,
           0.0,
           1.02166933,
           0.0,
           0.0,
           0.0,
           1.38132294,
           0.0,
           0.0,
           2.44530231,
           0.0}};
}