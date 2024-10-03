#include "common/execinfo.hpp"
#include <cassert>
#include <iostream>
#include <omp.h>
#include <rt_init.hpp>

static ExecInfo c_test_set_openmp_threads(int rank, int size, int n_threads)
{
  ExecInfo info;
  Core::SimulationParameters params;
  params.user_params.n_thread = n_threads;
  set_n_thread_current_rank(rank, size, info, params.user_params);
  std::cerr << "Number of threads: " << info.thread_per_process << '\n';
  return info;
}

// Test case 1: We have 4 nodes with 8 threads, so each node should have 2
// threads
void test_set_openmp_threads()
{
  int rank = 2;
  int size = 4;
  ExecInfo info = c_test_set_openmp_threads(rank, size, 8);
  assert(info.thread_per_process == 2); // Assertion
}

// Test case 2: One node with 16 threads, should have number of threads equal to
// the number of CPU cores
void test_set_openmp_threads_2()
{
  // int rank = 0;
  // int size = 1;
  // size_t n_core = omp_get_num_procs();
  // ExecInfo info = c_test_set_openmp_threads(rank, size, 16);
  // assert(info.thread_per_process == n_core); // Assertion
}

// Test case 3: 5 nodes and we are in the last (4), as 16%5 = 1, the last should
// have int(16/5 + 1) = 4 threads
void test_set_openmp_threads_3()
{
  int rank = 4;
  int size = 5;
  ExecInfo info = c_test_set_openmp_threads(rank, size, 16);
  assert(info.thread_per_process == 4); // Assertion
}

// Test case 4: Same situation as test 3 but we are not on the last node, so we
// have int(16/5) = 3 threads
void test_set_openmp_threads_4()
{
  int rank = 3;
  int size = 5;
  ExecInfo info = c_test_set_openmp_threads(rank, size, 16);
  assert(info.thread_per_process == 3); // Assertion
}

// Test case 5: Only one node with 4 threads, should have 4 threads per process
void test_set_openmp_threads_5()
{
  int rank = 0;
  int size = 1;
  ExecInfo info = c_test_set_openmp_threads(rank, size, 4);
  assert(info.thread_per_process == 4); // Assertion
}

// Test case 6: 2 nodes and we are in the first (0), so it should have int(8/2)
// = 4 threads
void test_set_openmp_threads_6()
{
  int rank = 0;
  int size = 2;
  ExecInfo info = c_test_set_openmp_threads(rank, size, 8);
  assert(info.thread_per_process == 4); // Assertion
}

// Test case 7: 3 nodes and we are in the second (1), so it should have
// int(12/3) = 4 threads
void test_set_openmp_threads_7()
{
  int rank = 1;
  int size = 3;
  ExecInfo info = c_test_set_openmp_threads(rank, size, 12);
  assert(info.thread_per_process == 4); // Assertion
}

int main()
{
  // Run the test
  test_set_openmp_threads();
  test_set_openmp_threads_2();
  test_set_openmp_threads_3();
  test_set_openmp_threads_4();
  test_set_openmp_threads_5();
  test_set_openmp_threads_6();
  test_set_openmp_threads_7();
  return 0;
}