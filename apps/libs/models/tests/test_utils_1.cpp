#include <models/utils.hpp>

#define ASSERT_EQUALS(a, b) \
    if ((a) != (b)) { \
        std::cerr << "Test failed: " << #a << " != " << #b << "\n"; \
    } else { \
        std::cout << "Test passed: " << #a << " == " << #b << "\n"; \
    }

int main() {
    // Test cases to validate min_var
    double result;

    // Test 1: Minimum of two numbers
    result = Models::min_var(3.0, 1.0);
    ASSERT_EQUALS(result, 1.0);

    // Test 2: Minimum of three numbers
    result = Models::min_var(7.5, 2.3, 9.8);
    ASSERT_EQUALS(result, 2.3);

    // Test 3: Minimum of an empty initializer list (corner case)
    result = Models::min_var();
    ASSERT_EQUALS(result, 0.);

    // Test 4: Minimum of a mix of positive and negative numbers
    result = Models::min_var(-5.0, 4.0, -3.0, 0.0);
    ASSERT_EQUALS(result, -5.0);

    // Test 5: Minimum of floating point numbers
    result = Models::min_var(0.5, 0.2, 0.8);
    ASSERT_EQUALS(result, 0.2);

    // Test 6: Minimum of integers
    result = Models::min_var(10, 20, 5, 15);
    ASSERT_EQUALS(result, 5);

    // Test 7: Minimum of large numbers
    result = Models::min_var(1000000, 999999, 1000001);
    ASSERT_EQUALS(result, 999999);

    return 0;
}
