#include <stdint.h>
#include <stdio.h>

#define NTRUP_P 761
#define NTRUP_Q 4591

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

void generate_omegas(int16_t *omega_powers, int16_t *omega_inv_powers, int16_t n, int16_t q) {
    int16_t primitive_root = 3; // Example primitive root, you may need to find the correct one
    int16_t omega = mod_exp(primitive_root, (q - 1) / n, q);
    int16_t omega_inv = mod_exp(omega, q - 2, q);

    omega_powers[0] = 1;
    omega_inv_powers[0] = 1;
    for (int i = 1; i < n; i++) {
        omega_powers[i] = (omega_powers[i - 1] * omega) % q;
        omega_inv_powers[i] = (omega_inv_powers[i - 1] * omega_inv) % q;
    }
}

int main() {
    int16_t omega_powers[NTRUP_P];
    int16_t omega_inv_powers[NTRUP_P];

    generate_omegas(omega_powers, omega_inv_powers, NTRUP_P, NTRUP_Q);

    FILE *file = fopen("omegas.txt", "w");
    if (file == NULL) {
        fprintf(stderr, "Error opening file for writing\n");
        return 1;
    }

    fprintf(file, "omega_powers[NTRUP_P] = { ");
    for (int i = 0; i < NTRUP_P; i++) {
        fprintf(file, "%d", omega_powers[i]);
        if (i < NTRUP_P - 1) {
            fprintf(file, ", ");
        }
    }
    fprintf(file, " };\n");

    fprintf(file, "omega_inv_powers[NTRUP_P] = { ");
    for (int i = 0; i < NTRUP_P; i++) {
        fprintf(file, "%d", omega_inv_powers[i]);
        if (i < NTRUP_P - 1) {
            fprintf(file, ", ");
        }
    }
    fprintf(file, " };\n");

    fclose(file);

    return 0;
}

