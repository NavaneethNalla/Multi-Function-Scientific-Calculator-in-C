#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <Windows.h>
typedef struct complex
{
    double real, imaginary;
} complex;

typedef struct
{
    float x, y, z;
} vector;
void show_animation(char str[]);
double area_circle(float radius);
double area_rectangle(float length, float breadth); // functions for area
double area_square(float side);
double area_triangle(float side1, float side2, float side3);

double perimeter_of_circle(float radius);
double perimeter_of_rectangle(float length, float breadth); // functions for perimeter
double perimeter_of_square(float side);
double perimeter_of_triangle(float side1, float side2, float side3);

long long decimal_to_binary(long long n); // functions for decimal and binary
long long binary_to_decimal(int n);

long long permutation(int n, int r); // permutation and combination functions
long long combination(int n, int r);

void vector_multiplication(vector a, vector b); // to calulate dot product and vector product

double logarithm_value(double x, double base); // logarithm function

#define max 10
// function to input a matrix
void input_of_matrices(int mat[max][max], int rows, int cols);
// function to display a matrix
void display_of_matrices(int mat[max][max], int rows, int cols);
// function to add n matrices
void addition_of_matrices(int mats[][max][max], int n, int result[max][max], int rows, int cols);
// function to multiply 2 matrices
void multiplication_of_2_matrices(int mat1[max][max], int mat2[max][max], int result[max][max], int rows1, int cols1, int rows2, int cols2);

void solution_of_quadratic_equation(float a, float b, float c); // quadratic equation solver

int is_prime(int n); // prime numbers check

// gcd of 2 numbers
int gcd_of_2_numbers(int a, int b);
// gcd of array of numbers
int gcd_of_n_numbers(int arr[], int n);
// lcm of two numbers
int lcm_of_2_numbers(int a, int b);
// lcm of array of numbers
int lcm_of_n_numbers(int arr[], int n);

int nth_term_in_fibonacci(int n);
long long factorial(int n);

void solution_of_linear_equations(double a1, double b1, double c1, double a2, double b2, double c2); // solution of linear equations

// differentiation of polynomial
void differentiate(int coeffs[], int degree);
// addition two polynomials
void sum_of_2_polynomials(int poly1[], int poly2[], int degree1, int degree2);
// integral of a polynomial
void integrate(int coeffs[], int degree);

// Function to print expansion of sin(x) up to n terms
void expansion_of_sine(int n);
// Function to print expansion of cos(x) up to n terms
void expansion_of_cos(int n);
// Function to print expansion of log(1+x) up to n terms
void expansion_of_log(int n);
// Function to print expansion of e^x up to n terms
void expansion_of_exponential(int n);

double mean_of_n_entries(double data[], int n);
double slope_of_linear_regression(double x[], double y[], int n); // linear regression
double intercept_of_linear_regression(double x[], double y[], int n, double slope);

void multiplication_of_2_complex_numbers(complex a, complex b, complex *result);
// multiplication of two complex numbers
void multiplication_of_n_complex_numbers(complex arr[], int n, complex *result);
// divide two complex numbers
void divison_of_2_complex_numbers(complex a, complex b, complex *result);

void bubbleSort(int array[], int n);

// Function to find the nth term of an arithmetic progression
float nth_term_of_ap(float first_term, float common_difference, int n);

// Function to find the sum of the first n terms of an arithmetic progression
float sum_of_nterms_of_ap(float first_term, float common_difference, int n);

// Function to check if a given term belongs to an arithmetic progression
int term_belongs_to_ap(float first_term, float common_difference, float term);

// function to print first 25 terms of an AP
void print_ap(float first_term, float common_difference);

// function to print first 25 terms of an GP
void print_gp(float first_term, float common_ratio);
// Function to find the nth term of an geometric progression
float nth_term_of_gp(float first_term, float common_ratio, int n);

// Function to find the sum of the first n terms of an geometric progression
float sum_of_nterms_of_gp(float first_term, float common_ratio, int n);

// Function to check if a given term belongs to an geometric progression
int term_belongs_to_gp(float first_term, float common_ratio, float term);

const double pi = 3.141592;
// Function to determine conversion factor
double conversion_length(char from[], char to[]);
double conversion_weight(char from[], char to[]);

// Function to calculate the mean of an array of numbers
float mean_of_n_numbers(int arr[], int n);

// Function to calculate the mode of an array of numbers
void calculate_mode(int arr[], int n);

// Function to calculate the variance of an array of numbers
float variance(int arr[], int n);

// Function to calculate the standard deviation of an array of numbers
float standard_deviation(int arr[], int n);

void bubbleSort(int array[], int n);
// Function to find determinant of a 2x2 matrix
double determinant(double a, double b, double c, double d);
int main()
{

        show_animation("\n\t\t\t--------------Welcome to calculator project in C-------------\n\n");
    int m = 1;
    do
    {
        system("cls");
        printf("1 for area                               2 for perimeter                  3 for binary to decimal conversion    4 for decimal to binary conversion\n");
        printf("5 for permutation                        6 for combination                7 for logarithm                       8 for matrices function\n9 for quadratic solver                   ");
        printf("10 for prime number              11 for basic calculations             12 for factorial\n13 for nth_term_in_fibonacci             14 for GCD                       15 for LCM                            ");
        printf("16 for conversion of units\n17 for polynomial differentiation        18 for polynomial integration    19 for polynomial adition             ");
        printf("20 for expansions\n21 for solution of 2 linear equatuons    22 for sorting                   23 for volume                         24 for linear regression\n25 for AP                                26 for GP                        27 for statistics                     ");
        printf("28 for vector multiplication\n");
        int choice;
        printf("enter your choice:\n");
        scanf("%d", &choice);
        switch (choice)
        {
        case 1:;
            int choice2;
            printf("1 for circle\n2 for rectangle\n3 for square\n4 for triangle\n");
            scanf("%d", &choice2);
            switch (choice2)
            {
            case 1:;
                float radius;
                printf("enter radius\n");
                scanf("%f", &radius);
                printf("area of circle is %lf\n", area_circle(radius));
                break;
            case 2:;
                float length, breadth;
                printf("enter length and breadth\n");
                scanf("%f%f", &length, &breadth);
                printf("area of rectangle is %lf\n", area_rectangle(length, breadth));
                break;
            case 3:;
                float side;
                printf("enter side\n");
                scanf("%f", &side);
                printf("the area of square is %lf\n", area_square(side));
                break;
            case 4:;
                float side1, side2, side3;
                printf("enter sides of triangle\n");
                scanf("%f%f%f", &side1, &side2, &side3);
                printf("the area of triangle is %lf\n", area_triangle(side1, side2, side3));
                break;
            }

            break;

        case 2:

            ;
            int choice3;
            printf("1 for circle\n2 for rectangle\n3 for square\n4 for triangle\n");
            scanf("%d", &choice3);
            switch (choice3)
            {
            case 1:;
                float radius;
                printf("enter radius\n");
                scanf("%f", &radius);
                printf("area of circle is %lf\n", perimeter_of_circle(radius));
                break;
            case 2:;
                float length, breadth;
                printf("enter length and breadth\n");
                scanf("%f%f", &length, &breadth);
                printf("area of rectangle is %lf\n", perimeter_of_rectangle(length, breadth));
                break;
            case 3:;
                float side;
                printf("enter side\n");
                scanf("%d", &side);
                printf("the area of square is %lf\n", perimeter_of_square(side));
                break;
            case 4:;
                float side1, side2, side3;
                printf("enter sides of triangle\n");
                scanf("%f%f%f", &side1, &side2, &side3);
                printf("the area of triangle is %lf\n", perimeter_of_triangle(side1, side2, side3));
                break;
            }
            break;
        case 3:;
            long long binary;
            printf("enter the number in binary\n");
            scanf("%lld", &binary);
            printf("the equivalent in decimal is %lld\n", binary_to_decimal(binary));
            break;
        case 4:;
            long long decimal;
            printf("enter the number in decimal\n");
            scanf("%lld", &decimal);
            printf("the equivalent in decimal is %lld\n", decimal_to_binary(decimal));
            break;
        case 5:;
            int n, r;
            printf("enter n and r for permutation(nPr)\n");
            scanf("%d%d", &n, &r);
            printf("the permutation is %lld\n", permutation(n, r));
            break;
        case 6:;
            int n1, r1;
            printf("enter n and r for combination(nCr)\n");
            scanf("%d%d", &n1, &r1);
            printf("the combination is %lld\n", combination(n1, r1));
            break;
        case 7:;
            double num, base;

            printf("Enter the number: ");
            scanf("%lf", &num);

            printf("Enter the base: ");
            scanf("%lf", &base);

            if (num <= 0 || base <= 0 || base == 1)
            {
                printf("Invalid input. Number and base must be positive and base must not be 1.\n");
            }
            else
            {
                double result = logarithm_value(num, base);
                printf("Logarithm of %.2lf to the base %.2lf is %.2lf\n", num, base, result);
            }
            break;
        case 8:;
            int choice4;
            printf("1 for addition\n2 for multiplication\n");
            scanf("%d", &choice4);
            switch (choice4)
            {
            case 1:
            {
                ;
                int n, rows, cols;
                printf("Enter the number of matrices: ");
                scanf("%d", &n);

                printf("Enter the number of rows and columns for the matrices: ");
                scanf("%d%d", &rows, &cols);

                int matrices[n][max][max];
                int result[max][max];

                // Input matrices
                for (int i = 0; i < n; i++)
                {
                    printf("Input for matrix %d:\n", i + 1);
                    input_of_matrices(matrices[i], rows, cols);
                }

                // Add matrices
                addition_of_matrices(matrices, n, result, rows, cols);

                // Display the result
                printf("Resultant matrix after addition:\n");
                display_of_matrices(result, rows, cols);
                break;
            }

            case 2:
            {
                ;

                int rows1, cols1, rows2, cols2;
                int mat1[max][max], mat2[max][max], result1[max][max];

                // Input dimensions and elements of the first matrix
                printf("Enter the number of rows and columns for the first matrix: ");
                scanf("%d%d", &rows1, &cols1);
                input_of_matrices(mat1, rows1, cols1);

                // Input dimensions and elements of the second matrix
                printf("Enter the number of rows and columns for the second matrix: ");
                scanf("%d%d", &rows2, &cols2);
                input_of_matrices(mat2, rows2, cols2);

                // Multiply matrices
                multiplication_of_2_matrices(mat1, mat2, result1, rows1, cols1, rows2, cols2);

                // Display the result
                printf("Resultant matrix after multiplication:\n");
                display_of_matrices(result1, rows1, cols2);
                break;
            }

            default:
                printf("Invalid choice. Please enter either 1 or 2.\n");
            }
            break;
        case 9:;
            float a, b, c;
            printf("enter a,b,c according to the equation ax^2+bx+c=0\n");
            scanf("%f%f%f", &a, &b, &c);
            solution_of_quadratic_equation(a, b, c);
            break;
        case 10:;
            int numn;
            printf("enter number\n");
            scanf("%d", &numn);
            (is_prime(numn) == 1) ? printf("the number is prime\n") : printf("the number is not prime\n");
            break;
        case 11:;
            int choice5;
            printf("1 for real numbers\n2 for complex numbers\n");
            scanf("%d", &choice5);
            switch (choice5)
            {
            case 1:
            {
                ;
                printf("1 for addition\n2 for multiplication\n3 for division\n4 for remainder\n");
                printf("enter your choice\n");
                int choice6;
                scanf("%d", &choice6);
                switch (choice6)
                {
                case 1:
                {
                    ;
                    int n2;
                    double sum = 0;
                    printf("enter number of entries\n");
                    scanf("%d", &n2);
                    double elements[n2];
                    for (int i = 1; i <= n2; i++)
                    {
                        printf("enter number:");
                        scanf("%lf", &elements[i - 1]);
                        sum += elements[i - 1];
                    }
                    printf("sum is %lf\n", sum);
                    break;
                }
                case 2:
                {
                    ;
                    int n3;
                    double multiplication = 1;
                    printf("enter number of entries\n");
                    scanf("%d", &n3);
                    double elements1[n3];
                    for (int i = 1; i <= n3; i++)
                    {
                        printf("enter number:\n");
                        scanf("%lf", &elements1[i - 1]);
                        multiplication *= elements1[i - 1];
                    }
                    printf("multiplication is %lf\n", multiplication);
                    break;
                }
                case 3:;
                    float a1, b1;
                    printf("enter the numerator and denominator\n");
                    scanf("%f%f", &a1, &b1);
                    float result = a1 / b1;
                    (a1 * b1 * result > 0) ? printf("result is %f\n", result) : printf("result is %f\n", -1 * result);
                    break;
                case 4:;
                    int a2;
                    int z;
                    printf("enter the numerator and denominator\n");
                    scanf("%d %d", &a2, &z);
                    int result2 = a2 % z;
                    printf("remainder is %f\n", result2);
                    break;
                }
                break;
            }
            case 2:
            {
                printf("1 for addition\n2 for multiplication\n3 for division\n");
                printf("enter your choice\n");
                int choice7;
                scanf("%d", &choice7);
                switch (choice7)
                {
                case 1:
                {
                    ;
                    int n4;
                    double a3 = 0, b3 = 0;
                    printf("enter number of entries:\n");
                    scanf("%d", &n4);
                    complex entry1[n4];
                    for (int i = 0; i < n4; i++)
                    {
                        printf("enter complex number %d\n:", i);
                        scanf("%lf %lf", &entry1[i].real, &entry1[i].imaginary);
                    }

                    for (int i = 0; i < n4; i++)
                    {
                        a3 += entry1[i].real;
                        b3 += entry1[i].imaginary;
                    }
                    printf("the sum is %lf+i%lf\n", a3, b3);
                    break;
                }
                case 2:
                {
                    ;
                    int n5;

                    printf("number of entries\n");
                    scanf("%d", &n5);
                    complex entries[n5];
                    printf("according to format x+iy , enter x and y\n");
                    for (int i = 0; i < n5; i++)
                    {
                        printf("enter entry number %d:\n", i + 1);
                        scanf("%lf%lf", &entries[i].real, &entries[i].imaginary);
                    }
                    complex result;
                    multiplication_of_n_complex_numbers(entries, n5, &result);
                    printf("the result of multiplication is %.2lf + %.2lfi\n", result.real, result.imaginary);
                    break;
                }
                case 3:
                {
                    ;
                    complex num1, num2, result;

                    // Input first complex number
                    printf("Enter the real part of the first complex number: ");
                    scanf("%lf", &num1.real);
                    printf("Enter the imaginary part of the first complex number: ");
                    scanf("%lf", &num1.imaginary);

                    // Input second complex number
                    printf("Enter the real part of the second complex number: ");
                    scanf("%lf", &num2.real);
                    printf("Enter the imaginary part of the second complex number: ");
                    scanf("%lf", &num2.imaginary);
                    divison_of_2_complex_numbers(num1, num2, &result);
                    // Display the result of division
                    printf("Result of division: %.2lf + %.2lfi\n", result.real, result.imaginary);
                    break;
                }
                }
                break;
            }
            }
            break;
        case 12:;
            int n6;
            printf("enter the number:\n");
            scanf("%d", &n6);
            printf("the factorial of %d is %d\n", n6, factorial(n6));
            break;
        case 13:;
            int n7;
            printf("enter the index of term in nth_term_in_fibonacci series\n");
            scanf("%d", &n7);
            printf("the %d term in nth_term_in_fibonacci series is %d\n", n7, nth_term_in_fibonacci(n7));
            break;
        case 14:
        {
            ;
            int numb_of_entries;
            printf("enter the number of entries:\n");
            scanf("%d", &numb_of_entries);
            int values[numb_of_entries];
            printf("enter the values:\n");
            for (int i = 0; i < numb_of_entries; i++)
            {
                scanf("%d", &values[i]);
            }
            int gcd = gcd_of_n_numbers(values, numb_of_entries);
            printf("GCD of the numbers is: %d\n", gcd);
            break;
        }
        case 15:
        {
            ;
            int num_of_entries;
            printf("enter the number of entries:\n");
            scanf("%d", &num_of_entries);
            int value[num_of_entries];
            printf("enter the values:\n");
            for (int i = 0; i < num_of_entries; i++)
            {
                scanf("%d", &value[i]);
            }
            int lcm = lcm_of_n_numbers(value, num_of_entries);
            printf("LCM of the numbers is: %d\n", lcm);
            break;
        }
        case 16:;
            int press;
            printf("1 for length conversion\n2 for weight conversion\n");
            scanf("%d", &press);
            switch (press)
            {
            case 1:
            {
                ;
                double value;
                char from[5], to[5];

                printf("Enter the value to convert:\n");
                scanf("%lf", &value);

                printf("Enter the unit to convert from (cm, m, km, dm, mm):\n");
                scanf("%s", from);

                printf("Enter the unit to convert to (cm, m, km, dm, mm):\n");
                scanf("%s", to);

                double factor = conversion_length(from, to);
                printf("%lf %s is equal to %lf %s.\n", value, from, value * factor, to);
                break;
            }
            case 2:
            {
                ;
                double value1;
                char from1[5], to1[5];

                printf("Enter the value to convert: ");
                scanf("%lf", &value1);

                printf("Enter the unit to convert from (g, kg, mg, lb, oz): ");
                scanf("%s", from1);

                printf("Enter the unit to convert to (g, kg, mg, lb, oz): ");
                scanf("%s", to1);

                double factor1 = conversion_weight(from1, to1);

                printf("%lf %s is equal to %lf %s.\n", value1, from1, value1 * factor1, to1);
                break;
            }
            }
            break;
        case 17:
        {
            ;
            int degree;
            printf("Enter the degree of the polynomial: ");
            scanf("%d", &degree);
            int coeffs[degree + 1];
            printf("Enter the coefficients of the polynomial (starting from the highest degree term): ");
            for (int i = 0; i <= degree; i++)
            {
                scanf("%d", &coeffs[i]);
            }
            differentiate(coeffs, degree);
            break;
        }
        case 18:
        {
            ;
            int degree;
            printf("Enter the degree of the polynomial: ");
            scanf("%d", &degree);
            int coeffs[degree + 1];
            printf("Enter the coefficients of the polynomial (starting from the constant term): ");
            for (int i = 0; i <= degree; i++)
            {
                scanf("%d", &coeffs[i]);
            }
            integrate(coeffs, degree);
            break;
        }
        case 19:
        {
            ;
            int degree1, degree2;
            printf("Enter the degree of the first polynomial: ");
            scanf("%d", &degree1);
            int poly1[degree1 + 1];
            printf("Enter the coefficients of the first polynomial (starting from the constant term): ");
            for (int i = 0; i <= degree1; i++)
            {
                scanf("%d", &poly1[i]);
            }
            printf("Enter the degree of the second polynomial: ");
            scanf("%d", &degree2);
            int poly2[degree2 + 1];
            printf("Enter the coefficients of the second polynomial (starting from the constant term): ");
            for (int i = 0; i <= degree2; i++)
            {
                scanf("%d", &poly2[i]);
            }
            sum_of_2_polynomials(poly1, poly2, degree1, degree2);
            break;
        }
        case 20:
        {
            ;
            int numb;
            char str1[] = "sinx";
            char str2[] = "cosx";
            char str3[] = "expo";
            char str4[] = "log";
            printf("enter number of terms to print in expansion\n");
            scanf("%d", &numb);
            char name[6];
            printf("enter expansion among sinx , cosx , expo for e^x ,log for log(1+x)\n");
            scanf("%s", name);
            if (!strcmp(name, str1))
            {
                expansion_of_sine(numb);
                printf("\n");
            }
            if (!strcmp(name, str2))
            {
                expansion_of_cos(numb);
                printf("\n");
            }
            if (!strcmp(name, str4))
            {
                expansion_of_log(numb);
                printf("\n");
            }
            if (!strcmp(name, str3))
            {
                expansion_of_exponential(numb);
                printf("\n");
            }
            if (c == 4)
            {
                printf("wrong choice\n");
            }
            break;
        }
        case 21:;
            double a1, b1, c1, a2, b2, c2;
            printf("Enter the coefficients of the first equation (ax + by = c):\n");
            printf("a1: ");
            scanf("%lf", &a1);
            printf("b1: ");
            scanf("%lf", &b1);
            printf("c1: ");
            scanf("%lf", &c1);
            printf("Enter the coefficients of the second equation (ax + by = c):\n");
            printf("a2: ");
            scanf("%lf", &a2);
            printf("b2: ");
            scanf("%lf", &b2);
            printf("c2: ");
            scanf("%lf", &c2);
            solution_of_linear_equations(a1, b1, c1, a2, b2, c2);
            break;
        case 22:
        {
            ;
            int n9;
            printf("Enter the number of elements: ");
            scanf("%d", &n9);
            int numberrs[n9];
            printf("Enter %d numbers:\n", n9);
            for (int i = 0; i < n9; ++i)
            {
                scanf("%d", &numberrs[i]);
            }
            printf("Array before sorting:\n");
            for (int i = 0; i < n9; ++i)
            {
                printf("%d ", numberrs[i]);
            }
            printf("\n");
            bubbleSort(numberrs, n9);
            printf("Array after sorting using Bubble Sort:\n");
            for (int i = 0; i < n9; ++i)
            {
                printf("%d ", numberrs[i]);
            }
            printf("\n");
            break;
        }
        case 23:;
            int cchoice;
            printf("1 for cuboid\n2 for triangular prism\n3 for cylinder\n4 for pyramid\n5 for tetrahedron\n");
            scanf("%d", &cchoice);
            switch (cchoice)
            {
            case 1:;
                double length, width, heigght;
                printf("Enter the length, width, and height of the rectangular prism: ");
                scanf("%lf %lf %lf", &length, &width, &heigght);
                double resultt = 2 * (length * width + width * heigght + heigght * length);
                printf("Surface area of the rectangular prism: %.2f\n", resultt);
                break;
            case 2:;
                double base, heighht, prismHeight;
                printf("Enter the base, height, and prism height of the triangular prism: ");
                scanf("%lf %lf %lf", &base, &heighht, &prismHeight);
                double resullt = base * heighht * prismHeight / 2;
                printf("Volume of the triangular prism: %.2f\n", resullt);
                break;
            case 3:;
                double raddius, heigt, ressult;
                printf("enter the radius and height\n");
                scanf("%lf %lf", &raddius, &heigt);
                ressult = pi * pow(raddius, 2) * heigt;
                printf("the volume of the desired cylinder is %lf", ressult);
                break;
            case 4:;
                double basee, height;
                printf("Enter the base length and height of the pyramid: ");
                scanf("%lf %lf", &basee, &height);
                double reesult = (base * basee * height) / 3;
                printf("Volume of the pyramid: %.2f\n", reesult);
                break;
            case 5:;
                double side;
                printf("Enter the side length of the tetrahedron: ");
                scanf("%lf", &side);
                double resuult = (sqrt(3) * side * side * side) / 12;
                printf("Volume of the tetrahedron: %.2f\n", resuult);
                break;
            }
            break;
        case 24:;
            int n10;
            printf("Enter the number of data points: ");
            scanf("%d", &n10);

            double *x = (double *)malloc(n10 * sizeof(double));
            double *y = (double *)malloc(n10 * sizeof(double));

            printf("Enter the data points (x, y):\n");
            for (int i = 0; i < n10; i++)
            {
                scanf("%lf %lf", &x[i], &y[i]);
            }
            double slope = slope_of_linear_regression(x, y, n10);
            double intercept = intercept_of_linear_regression(x, y, n10, slope);
            printf("Regression Line Equation: y = %.2lfx + %.2lf\n", slope, intercept);
            free(x);
            free(y);
            break;
        case 25:;
            float first_term, common_difference, term;
            int nth;
            printf("Enter the first term of the AP: ");
            scanf("%f", &first_term);
            printf("Enter the common difference of the AP: ");
            scanf("%f", &common_difference);
            printf("Enter the value of n: ");
            scanf("%d", &nth);
            printf("Enter the term to check if it belongs to the AP: ");
            scanf("%f", &term);
            printf("The %dth term of the AP is: %.2f\n", nth, nth_term_of_ap(first_term, common_difference, nth));
            printf("The sum of the first %d terms of the AP is: %.2f\n", nth, sum_of_nterms_of_ap(first_term, common_difference, nth));
            if (term_belongs_to_ap(first_term, common_difference, term))
                printf("The term %.2f belongs to the AP.\n", term);
            else
                printf("The term %.2f does not belong to the AP.\n", term);
            print_ap(first_term, common_difference);

            break;
        case 26:

            ;
            float first_term1, common_ratio, term1;
            int nth2;
            printf("Enter the first term of the GP: ");
            scanf("%f", &first_term1);
            printf("Enter the common ratio of the AP: ");
            scanf("%f", &common_ratio);
            printf("Enter the value of n: ");
            scanf("%d", &nth2);
            printf("Enter the term to check if it belongs to the GP: ");
            scanf("%f", &term1);
            printf("The %dth term of the GP is: %.2f\n", nth2, nth_term_of_gp(first_term1, common_ratio, nth2));
            printf("The sum of the first %d terms of the AP is: %.2f\n", nth2, sum_of_nterms_of_gp(first_term1, common_ratio, nth2));
            if (term_belongs_to_gp(first_term1, common_ratio, term1))
                printf("The term %.2f belongs to the AP.\n", term1);
            else
                printf("The term %.2f does not belong to the AP.\n", term1);
            print_gp(first_term1, common_ratio);

            break;

        case 27:
        {
            ;
            int t;
            printf("Enter the number of elements: ");
            scanf("%d", &t);
            int *arrr = (int *)malloc(t * sizeof(int));

            printf("Enter the elements: ");
            for (int i = 0; i < t; i++)
            {
                scanf("%d", &arrr[i]);
            }

            float mean = mean_of_n_numbers(arrr, t);
            printf("Mean: %.2f\n", mean);

            printf("Mode: ");
            calculate_mode(arrr, t);

            float variaance = variance(arrr, t);
            printf("Variance: %.2f\n", variaance);

            float standard_deviiation = standard_deviation(arrr, t);
            printf("Standard Deviation: %.2f\n", standard_deviiation);
            break;
        }
        case 28:
        {
            ;
            vector A, B;
            printf("enter 1st vector as:\n");
            scanf("%f %f %f", &A.x, &A.y, &A.z);
            printf("enter 2nd vector as:\n");
            scanf("%f %f %f", &B.x, &B.y, &B.z);
            vector_multiplication(A, B);
            break;
        }

        default:
        {
            ;
            printf("wrong choice\n");
            break;
        }
        }

        printf("press 1 for another operation and 0 to exit\n");
        scanf("%d", &m);
    } while (m != 0);
    return 0;
}
void multiplication_of_2_complex_numbers(complex a, complex b, complex *result)
{
    result->real = (a.real * b.real) - (a.imaginary * b.imaginary);
    result->imaginary = (a.real * b.imaginary) + (a.imaginary * b.real);
}
void multiplication_of_n_complex_numbers(complex arr[], int n, complex *result)
{
    result->real = 1;
    result->imaginary = 0;
    for (int i = 0; i < n; i++)
    {
        complex temp;
        multiplication_of_2_complex_numbers(*result, arr[i], &temp);
        *result = temp;
    }
}
// Function to divide two complex numbers
void divison_of_2_complex_numbers(complex a, complex b, complex *result)
{
    double divisor = b.real * b.real + b.imaginary * b.imaginary;
    result->real = (a.real * b.real + a.imaginary * b.imaginary) / divisor;
    result->imaginary = (a.imaginary * b.real - a.real * b.imaginary) / divisor;
}
// Function to determine conversion factor
double conversion_length(char from[], char to[])
{
    // Conversion factors
    if (strcmp(from, "m") == 0 && strcmp(to, "cm") == 0)
        return 100; // 1 meter = 100 centimeters
    else if (strcmp(from, "cm") == 0 && strcmp(to, "m") == 0)
        return 0.01; // 1 centimeter = 0.01 meters
    else if (strcmp(from, "m") == 0 && strcmp(to, "km") == 0)
        return 0.001; // 1 meter = 0.001 kilometers
    else if (strcmp(from, "km") == 0 && strcmp(to, "m") == 0)
        return 1000; // 1 kilometer = 1000 meters
    else if (strcmp(from, "cm") == 0 && strcmp(to, "km") == 0)
        return 0.00001; // 1 centimeter = 0.00001 kilometers
    else if (strcmp(from, "km") == 0 && strcmp(to, "cm") == 0)
        return 100000; // 1 kilometer = 100000 centimeters
    else if (strcmp(from, "m") == 0 && strcmp(to, "mm") == 0)
        return 1000; // 1 meter = 1000 millimeters
    else if (strcmp(from, "mm") == 0 && strcmp(to, "m") == 0)
        return 0.001; // 1 millimeter = 0.001 meters
    else if (strcmp(from, "cm") == 0 && strcmp(to, "mm") == 0)
        return 10; // 1 centimeter = 10 millimeters
    else if (strcmp(from, "mm") == 0 && strcmp(to, "cm") == 0)
        return 0.1; // 1 millimeter = 0.1 centimeters
    else if (strcmp(from, "dm") == 0 && strcmp(to, "mm") == 0)
        return 100; // 1 decimeter = 100 millimeters
    else if (strcmp(from, "mm") == 0 && strcmp(to, "dm") == 0)
        return 0.01; // 1 millimeter = 0.01 decimeters
    else if (strcmp(from, "m") == 0 && strcmp(to, "dm") == 0)
        return 10; // 1 meter = 10 decimeters
    else if (strcmp(from, "cm") == 0 && strcmp(to, "dm") == 0)
        return 0.1; // 1 centimeter = 0.1 decimeters
    else if (strcmp(from, "km") == 0 && strcmp(to, "dm") == 0)
        return 10000; // 1 kilometer = 10000 decimeters
    else if (strcmp(from, "km") == 0 && strcmp(to, "mm") == 0)
        return 1000000; // 1 kilometer = 1000000 milimeters
    else if (strcmp(from, "mm") == 0 && strcmp(to, "km") == 0)
        return 0.000001; // 1 millimeter = 0.000001 kilometers
    else if (strcmp(from, "dm") == 0 && strcmp(to, "m") == 0)
        return 0.1; // 1 decimeter = 0.01 meters
    else if (strcmp(from, "dm") == 0 && strcmp(to, "km") == 0)
        return 0.0001; // 1 decimeter = 0.0001 kilometers
    else if (strcmp(from, "dm") == 0 && strcmp(to, "cm") == 0)
        return 10; // 1 decimeter = 10 centimeters
    else
        return 1; // No conversion needed
}
// Function to calculate the derivative of a polynomial
void differentiate(int coeffs[], int degree)
{
    printf("The derivative of the polynomial is: ");
    for (int i = 0; i < degree; i++)
    {
        printf("%dx^%d ", coeffs[i] * (degree - i), degree - i - 1);
        if (i < degree - 1)
            printf("+ ");
    }
    printf("\n");
}
// Function to add two polynomials
void sum_of_2_polynomials(int poly1[], int poly2[], int degree1, int degree2)
{
    int maxDegree = (degree1 > degree2) ? degree1 : degree2;
    int sum[maxDegree + 1];
    for (int i = 0; i <= maxDegree; i++)
    {
        sum[i] = 0;
    }
    for (int i = 0; i <= degree1; i++)
    {
        sum[i] += poly1[i];
    }
    for (int i = 0; i <= degree2; i++)
    {
        sum[i] += poly2[i];
    }
    printf("The sum of the polynomials is: ");
    for (int i = maxDegree; i >= 0; i--)
    {
        if (sum[i] != 0)
        {
            printf("%dx^%d ", sum[i], i);
            if (i > 0)
                printf("+ ");
        }
    }
    printf("\n");
}
// Function to calculate the integral of a polynomial
void integrate(int coeffs[], int degree)
{
    printf("The integral of the polynomial is: ");
    for (int i = 0; i <= degree; i++)
    {
        if (coeffs[i] != 0)
        {
            printf("%d/%d*x^%d ", coeffs[i], i + 1, i + 1);
            if (i < degree)
                printf("+ ");
        }
    }
    printf("+ C\n"); // Add constant of integration
}

int nth_term_in_fibonacci(int n)
{
    if (n <= 2)
        return n - 1;

    return nth_term_in_fibonacci(n - 1) + nth_term_in_fibonacci(n - 2);
}

long long factorial(int n)
{
    long long ans = 1;
    for (int i = 1; i <= n; i++)
    {
        ans *= i;
    }
    return ans;
}

double area_circle(float radius)
{
    double area = pi * pow(radius, 2);
    return area;
}
double area_rectangle(float length, float breadth)
{
    double area = length * breadth;
    return area;
}
double area_square(float side)
{
    return (double)(pow(side, 2));
}
double area_triangle(float side1, float side2, float side3)
{
    double s = (side1 + side2 + side3) / 2;
    double area_square = s * (s - side1) * (s - side2) * (s - side3);
    return pow(area_square, 0.5);
}
double perimeter_of_circle(float radius)
{
    return (2 * pi * radius);
}
double perimeter_of_rectangle(float length, float breadth)
{
    return 2 * (length + breadth);
}
double perimeter_of_square(float side)
{
    return 4 * side;
}
double perimeter_of_triangle(float side1, float side2, float side3)
{
    return (side1 + side2 + side3);
}

long long decimal_to_binary(long long n)
{
    int x = 1;
    int ans = 0;
    while (x <= n)
    {
        x *= 2;
    }
    while (x > 0)
    {
        int lastdigit = n / x;
        n -= lastdigit * x;
        ans = ans * 10 + lastdigit;
        x /= 2;
    }
    return ans;
}
long long binary_to_decimal(int n)
{
    int x = 1;
    int ans = 0;

    while (n > 0)
    {
        int y = n % 10;
        ans += x * y;
        x *= 2;
        n = n / 10;
    }
    return ans;
}
long long permutation(int n, int r)
{
    int ans = factorial(n) / factorial(n - r);
    return ans;
}
long long combination(int n, int r)
{
    int ans = factorial(n) / (factorial(r) * factorial(n - r));
    return ans;
}
// Function to find GCD of two numbers
int gcd_of_2_numbers(int a, int b)
{
    if (b == 0)
        return a;
    return gcd_of_2_numbers(b, a % b);
}

// Function to find GCD of array of numbers
int gcd_of_n_numbers(int arr[], int n)
{
    int result = arr[0];
    for (int i = 1; i < n; i++)
    {
        result = gcd_of_2_numbers(result, arr[i]);
    }
    return result;
}

// Function to find LCM of two numbers
int lcm_of_2_numbers(int a, int b)
{
    return (a * b) / gcd_of_2_numbers(a, b);
}

// Function to find LCM of array of numbers
int lcm_of_n_numbers(int arr[], int n)
{
    int result = arr[0];
    for (int i = 1; i < n; i++)
    {
        result = lcm_of_2_numbers(result, arr[i]);
    }
    return result;
}

double logarithm_value(double x, double base)
{
    return log(x) / log(base);
}
void input_of_matrices(int mat[max][max], int rows, int cols)
{
    printf("Enter the elements of the matrix (%dx%d):\n", rows, cols);
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            printf("Enter element [%d][%d]: ", i + 1, j + 1);
            scanf("%d", &mat[i][j]);
        }
    }
}

// Function to display a matrix
void display_of_matrices(int mat[max][max], int rows, int cols)
{
    printf("The matrix is:\n");
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            printf("%d ", mat[i][j]);
        }
        printf("\n");
    }
}

// Function to add n matrices
void addition_of_matrices(int mats[][max][max], int n, int result[max][max], int rows, int cols)
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            result[i][j] = 0;
            for (int k = 0; k < n; k++)
            {
                result[i][j] += mats[k][i][j];
            }
        }
    }
}

void multiplication_of_2_matrices(int mat1[max][max], int mat2[max][max], int result[max][max], int rows1, int cols1, int rows2, int cols2)
{
    if (cols1 != rows2)
    {
        printf("Matrix multiplication not possible: Number of columns in the first matrix "
               "is not equal to the number of rows in the second matrix.\n");
        return;
    }

    for (int i = 0; i < rows1; i++)
    {
        for (int j = 0; j < cols2; j++)
        {
            result[i][j] = 0;
            for (int k = 0; k < cols1; k++)
            {
                result[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
    }
}
void solution_of_quadratic_equation(float a, float b, float c)
{
    int D = b * b - 4 * a * c;
    float root1, root2, x, y;
    if (D > 0)
    {
        root1 = ((-b) + pow(D, 0.5)) / (2 * a);
        root2 = ((-b) - pow(D, 0.5)) / (2 * a);
    }
    else
    {
        x = (-b) / (2 * a);
        y = pow(-D, 0.5) / (2 * a);
    }
    (D > 0) ? printf("the roots of the equation are %f and %f\n", root1, root2) : printf("the roots of the equation are %f+i%f and %f-i%f\n", x, y);
}
int is_prime(int n)
{

    if (n == 2)
    {
        return 1;
    }
    int count = 0;
    for (int i = 2; i < n; i++)
    {
        if (n % i == 0)
        {
            count++;
            break;
        }
    }
    return (count == 0) ? 1 : 0;
}
// Function to print expansion of sin(x) up to n terms
void expansion_of_sine(int n)
{
    printf("Expansion of sin(x) up to %d terms:\n", n);
    for (int i = 0; i < n; ++i)
    {
        printf("%.2f*x^%d/%d! ", pow(-1, i), 2 * i + 1, 2 * i + 1);
        if (i < n - 1)
            printf("+ ");
    }
    printf("\n");
}

// Function to print expansion of cos(x) up to n terms
void expansion_of_cos(int n)
{
    printf("Expansion of cos(x) up to %d terms:\n", n);
    for (int i = 0; i < n; ++i)
    {
        printf("%.2f*x^%d/%d! ", pow(-1, i), 2 * i, 2 * i);
        if (i < n - 1)
            printf("+ ");
    }
    printf("\n");
}

// Function to print expansion of log(1+x) up to n terms
void expansion_of_log(int n)
{
    printf("Expansion of log(1+x) up to %d terms:\n", n);
    for (int i = 1; i <= n; ++i)
    {
        printf("%.2f*x^%d/%d ", pow(-1, i + 1), i, i);
        if (i < n)
            printf("+ ");
    }
    printf("\n");
}

// Function to print expansion of e^x up to n terms
void expansion_of_exponential(int n)
{
    printf("Expansion of e^x up to %d terms:\n", n);
    for (int i = 0; i < n; ++i)
    {
        printf("x^%d/%d! ", i, i);
        if (i < n - 1)
            printf("+ ");
    }
    printf("\n");
}
// Function to find determinant of a 2x2 matrix
double determinant(double a, double b, double c, double d)
{
    return a * d - b * c;
}

// Function to solve a system of linear equations
void solution_of_linear_equations(double a1, double b1, double c1, double a2, double b2, double c2)
{
    double detA = determinant(a1, b1, a2, b2);
    double detX = determinant(b1, c1, b2, c2);
    double detY = determinant(c1, a1, c2, a2);

    if (detA != 0)
    {
        double x = detX / detA;
        double y = detY / detA;
        printf("Solution: x = %.2f, y = %.2f\n", x, y);
    }
    else
    {
        if (detX == 0 && detY == 0)
            printf("Infinite solutions.\n");
        else
            printf("No solution.\n");
    }
}
void bubbleSort(int array[], int n)
{
    int i, j;
    for (i = 0; i < n - 1; i++)
    {
        for (j = 0; j < n - i - 1; j++)
        {
            if (array[j] > array[j + 1])
            {
                // Swap array[j] and array[j+1]
                int temp = array[j];
                array[j] = array[j + 1];
                array[j + 1] = temp;
            }
        }
    }
}
// Function to calculate the mean_of_n_entries of an array
double mean_of_n_entries(double data[], int n)
{
    double sum = 0.0;
    for (int i = 0; i < n; i++)
    {
        sum += data[i];
    }
    return sum / n;
}

// Function to calculate the slope (m) of the regression line
double slope_of_linear_regression(double x[], double y[], int n)
{
    double sum_xy = 0.0, sum_x = 0.0, sum_y = 0.0, sum_x_squared = 0.0;
    for (int i = 0; i < n; i++)
    {
        sum_xy += x[i] * y[i];
        sum_x += x[i];
        sum_y += y[i];
        sum_x_squared += x[i] * x[i];
    }
    return (n * sum_xy - sum_x * sum_y) / (n * sum_x_squared - sum_x * sum_x);
}

// Function to calculate the intercept (b) of the regression line
double intercept_of_linear_regression(double x[], double y[], int n, double slope)
{
    double x_mean_of_n_entries = mean_of_n_entries(x, n);
    double y_mean_of_n_entries = mean_of_n_entries(y, n);
    return y_mean_of_n_entries - slope * x_mean_of_n_entries;
}
// Function to find the nth term of an arithmetic progression
float nth_term_of_ap(float first_term, float common_difference, int n)
{
    return first_term + (n - 1) * common_difference;
}

// Function to find the sum of the first n terms of an arithmetic progression
float sum_of_nterms_of_ap(float first_term, float common_difference, int n)
{
    return (n / 2.0) * (2 * first_term + (n - 1) * common_difference);
}

// Function to check if a given term belongs to an arithmetic progression
int term_belongs_to_ap(float first_term, float common_difference, float term)
{
    if (((term - first_term) / common_difference) == ((int)((term - first_term) / common_difference)))
        return 1; // Term belongs to the AP
    else
        return 0; // Term does not belong to the AP
}
// Function to find the nth term of an geometric progression
float nth_term_of_gp(float first_term, float common_ratio, int n)
{
    return first_term * pow(common_ratio, n - 1);
}

// Function to find the sum of the first n terms of an geometric progression
float sum_of_nterms_of_gp(float first_term, float common_ratio, int n)
{
    return first_term * ((1 - pow(common_ratio, n)) / (1 - common_ratio));
}

// Function to check if a given term belongs to an geometric progression
int term_belongs_to_gp(float first_term, float common_ratio, float term)
{
    if (logarithm_value(term / first_term, common_ratio) == (int)logarithm_value(term / first_term, common_ratio))
    {
        return 1;
    } // Term belongs to the AP
    else
    {
        return 0;
    } // Term does not belong to the AP
}
void print_ap(float first_term, float common_difference)
{
    for (int i = 1; i < 25; i++)
    {
        printf("%f  ", first_term + (i - 1) * common_difference);
    }
    printf(". . . .\n");
}
void print_gp(float first_term, float common_ratio)
{
    for (int i = 1; i < 25; i++)
    {
        printf("%f  ", first_term * pow(common_ratio, i - 1));
    }
    printf(". . . .\n");
}

float mean_of_n_numbers(int arr[], int n)
{
    float sum = 0;
    for (int i = 0; i < n; i++)
    {
        sum += arr[i];
    }
    return sum / n;
}

void calculate_mode(int arr[], int n)
{
    int b = 0, number_of_modes = 0;
    int mode[n];
    for (int i = 0; i < n; i++)
    {
        int count = 1;
        for (int j = i + 1; j < n; j++)
        {
            if (arr[i] == arr[j])
            {
                count++;
            }
        }
        if (count > b)
        {
            b = count;
            number_of_modes = 1;
            mode[0] = arr[i];
        }
        else if (count == b)
        {
            mode[number_of_modes] = arr[i];
            number_of_modes++;
        }
    }
    if (b == 1)
    {
        printf("No mode\n");
    }
    else
    {
        printf("Mode(s): ");
        for (int i = 0; i < number_of_modes; i++)
        {
            printf("%d ", mode[i]);
        }
        printf("\n");
    }
}
float variance(int arr[], int n)
{
    float mean = mean_of_n_numbers(arr, n);
    float variance = 0;
    for (int i = 0; i < n; i++)
    {
        variance += pow(arr[i] - mean, 2);
    }
    return variance / n;
}
float standard_deviation(int arr[], int n)
{
    return sqrt(variance(arr, n));
}
void vector_multiplication(vector a, vector b)
{
    float dot_product = (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
    vector cross_product;
    cross_product.x = a.y * b.z - a.z * b.y;
    cross_product.y = a.z * b.x - a.x * b.z;
    cross_product.z = a.x * b.y - a.y * b.x;
    printf("the dot product of (%f,%f,%f)  and  (%f,%f,%f) is %f\n", a.x, a.y, a.z, b.x, b.y, b.z, dot_product);
    printf("the cross product of (%f,%f,%f)  and  (%f,%f,%f) is  (%f,%f,%f) \n", a.x, a.y, a.z, b.x, b.y, b.z, cross_product.x, cross_product.y, cross_product.z);
}
void show_animation(char str[])
{
    for (int i = 0; i < strlen(str); i++)
    {
        Sleep(50);
        printf("%c", str[i]);
    }
    return;
}
double conversion_weight(char from[], char to[])
{
    // Conversion factors
    if (strcmp(from, "g") == 0 && strcmp(to, "kg") == 0)
        return 0.001;
    else if (strcmp(from, "g") == 0 && strcmp(to, "mg") == 0)
        return 1000;
    else if (strcmp(from, "g") == 0 && strcmp(to, "lb") == 0)
        return 0.0022;
    else if (strcmp(from, "g") == 0 && strcmp(to, "oz") == 0)
        return 0.353;
    else if (strcmp(from, "kg") == 0 && strcmp(to, "g") == 0)
        return 1000;
    else if (strcmp(from, "kg") == 0 && strcmp(to, "mg") == 0)
        return 1000000;
    else if (strcmp(from, "kg") == 0 && strcmp(to, "lb") == 0)
        return 2.205;
    else if (strcmp(from, "kg") == 0 && strcmp(to, "oz") == 0)
        return 35.274;
    else if (strcmp(from, "mg") == 0 && strcmp(to, "g") == 0)
        return .001;
    else if (strcmp(from, "mg") == 0 && strcmp(to, "kg") == 0)
        return 0.000001;
    else if (strcmp(from, "mg") == 0 && strcmp(to, "lb") == 0)
        return 0.000002205;
    else if (strcmp(from, "mg") == 0 && strcmp(to, "oz") == 0)
        return 0.00003527;
    else if (strcmp(from, "lb") == 0 && strcmp(to, "g") == 0)
        return 453.5;
    else if (strcmp(from, "lb") == 0 && strcmp(to, "kg") == 0)
        return 0.4544;
    else if (strcmp(from, "lb") == 0 && strcmp(to, "mg") == 0)
        return 453592.37;
    else if (strcmp(from, "lb") == 0 && strcmp(to, "oz") == 0)
        return 16;
    else if (strcmp(from, "oz") == 0 && strcmp(to, "g") == 0)
        return 28.35;
    else if (strcmp(from, "oz") == 0 && strcmp(to, "kg") == 0)
        return 0.0283;
    else if (strcmp(from, "oz") == 0 && strcmp(to, "lb") == 0)
        return 0.0625;
    else if (strcmp(from, "oz") == 0 && strcmp(to, "mg") == 0)
        return 28349.523;
    else
        return 1; // No conversion needed
}
