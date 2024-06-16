#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

// Function to perform modular exponentiation
int16_t mod_exp(int16_t base, int16_t exp, int16_t mod) {
    int16_t result = 1;
    while (exp > 0) {
        if (exp % 2 == 1) {
            result = (result * base) % mod;
        }
        base = (base * base) % mod;
        exp = exp / 2;
    }
    return result;
}

// Function to check if a candidate is a primitive root of unity
int is_primitive_root_of_unity(int16_t candidate, int16_t n, int16_t q) {
    if (mod_exp(candidate, n, q) != 1) {
        return 0;
    }
    for (int i = 1; i < n; i++) {
        if (n % i == 0 && mod_exp(candidate, n / i, q) == 1) {
            return 0;
        }
    }
    return 1;
}

// Function to find a primitive root of unity
int16_t find_primitive_root_of_unity(int16_t n, int16_t q) {
    for (int16_t candidate = 2; candidate < q; candidate++) {
        if (is_primitive_root_of_unity(candidate, n, q)) {
            return candidate;
        }
    }
    return -1;
}

// Function to compute the NTT using iterative method
void ntt(int16_t *a, int16_t n, int16_t omega, int16_t q) {
    for (int len = 1; len < n; len *= 2) {
        int16_t wlen = mod_exp(omega, n / (2 * len), q);
        for (int i = 0; i < n; i += 2 * len) {
            int16_t w = 1;
            for (int j = 0; j < len; j++) {
                int16_t u = a[i + j];
                int16_t v = (a[i + j + len] * w) % q;
                a[i + j] = (u + v) % q;
                a[i + j + len] = (u - v + q) % q;
                w = (w * wlen) % q;
            }
        }
    }
}

// Function to compute the inverse NTT using iterative method
void intt(int16_t *a, int16_t n, int16_t omega_inv, int16_t q) {
    ntt(a, n, omega_inv, q);
    int16_t inv_n = mod_exp(n, q - 2, q);
    for (int i = 0; i < n; i++) {
        a[i] = (a[i] * inv_n) % q;
    }
}

// Function to compute convolution using NTT
void convolution_ntt(int16_t *arr1, int16_t *arr2, int16_t q_ntt, int16_t q_target, int16_t n) {
    int16_t omega = find_primitive_root_of_unity(n, q_ntt);
    if (omega == -1) {
        printf("No suitable root of unity found\n");
        return;
    }
    int16_t omega_inv = mod_exp(omega, q_ntt - 2, q_ntt);

    // Pad arrays to length n
    int16_t *a1 = calloc(n, sizeof(int16_t));
    int16_t *a2 = calloc(n, sizeof(int16_t));
    memcpy(a1, arr1, 761 * sizeof(int16_t));
    memcpy(a2, arr2, 761 * sizeof(int16_t));

    // Compute NTT of both arrays
    ntt(a1, n, omega, q_ntt);
    ntt(a2, n, omega, q_ntt);

    // Pointwise multiplication
    int16_t *result_ntt = malloc(n * sizeof(int16_t));
    for (int i = 0; i < n; i++) {
        result_ntt[i] = (a1[i] * a2[i]) % q_ntt;
    }

    // Compute inverse NTT
    intt(result_ntt, n, omega_inv, q_ntt);

    // Reduce result modulo q_target
    int16_t *result = malloc(n * sizeof(int16_t));
    for (int i = 0; i < n; i++) {
        result[i] = result_ntt[i] % q_target;
    }

    // Print the result
    for (int i = 0; i < n; i++) {
        printf("%d ", result[i]);
    }
    printf("\n");

    // Free allocated memory
    free(a1);
    free(a2);
    free(result_ntt);
    free(result);
}

int main() {
    int16_t q_ntt = 674182657;
    int16_t q_target = 4591;
    int16_t n = 761;
    int16_t arr1[761];
    int16_t arr2[761];

    // Initialize arrays with example values
    for (int i = 0; i < 761; i++) {
        arr1[i] = 1;
        arr2[i] = 1;
    }

    convolution_ntt(arr1, arr2, q_ntt, q_target, n);
    return 0;
}

