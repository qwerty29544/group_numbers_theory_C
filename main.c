#include<stdio.h>
#include<malloc.h>
#include<math.h>


int prod_array_i(int *array, int length) {
    int result = 1;
    int i;
    for (i = 0; i < length; i++) {
        result *= array[i];
    }
    return result;
}


int sum_array_i(int *array, int length) {
    int result = 0;
    int i;
    for (i = 0; i < length; i++) {
        result += array[i];
    }
    return result;
}


short array_cmp_i(int *array1, int *array2, int length) {
    short result = 1;
    int i;
    for (i = 0; i < length; i++) {
        if (array1[i] != array2[i])
            result = 0;
    }
    return result;
}


typedef struct china_resid {
    int problem_length;
    int *residuals;
    int *modules_system;
    int modules_prod;
} crt;


crt init_crt(int problem_length,
             int *residuals,
             int *modules_system,
             int modules_prod) {
    crt crp_object;
    crp_object.problem_length = problem_length;
    crp_object.residuals = residuals;
    crp_object.modules_system = modules_system;
    if (modules_prod == 0) {
        crp_object.modules_prod = prod_array_i(crp_object.modules_system, crp_object.problem_length);
    } else {
        crp_object.modules_prod = modules_prod;
    }
    return crp_object;
};


crt init_crt_from_int(int problem_length,
                      int x,
                      int *modules_system,
                      int modules_prod) {
    crt crp_object;
    int *residuals_x = (int *)malloc(problem_length * sizeof(int));
    
    int i;
    for (i = 0; i < problem_length; i++) {
        residuals_x[i] = x % modules_system[i];
    }
    
    crp_object = init_crt(problem_length, residuals_x, modules_system, modules_prod);
    return crp_object;
}


void print_crt_problem(crt problem) {
    printf("\nChina residuals theorem problem:\n");
    
    int iter_index;
    for (iter_index = 0; 
         iter_index < problem.problem_length; 
         iter_index++) {
         printf("x = %4d{%d}\n", 
                problem.residuals[iter_index], 
                problem.modules_system[iter_index]);
    }
    
    printf("Find x \n");
};


void print_crt_problem_latex(crt problem) {
    printf("\n$$\n\\left\\{\n\\begin{matrix} \n");
    int iter_index;
    for (iter_index = 0; 
         iter_index < problem.problem_length - 1; 
         iter_index++) {
         printf("x = %4d\\{%d\\}\\\\\n", 
                problem.residuals[iter_index], 
                problem.modules_system[iter_index]);
    }
    printf("x = %4d\\{%d\\}\n", 
            problem.residuals[problem.problem_length - 1], 
            problem.modules_system[problem.problem_length - 1]);

    printf("\\end{matrix}\n\\right.\n$$\n");
}


void print_crt_number(crt number) {
    printf("Number vector = (%d", number.residuals[0]);
    
    int i;
    for (i = 1; i < number.problem_length; i++) {
        printf(", %d", number.residuals[i]);
    }
    
    printf("), modules vector = (%d", number.modules_system[0]);
    for (i = 1; i < number.problem_length; i++) {
        printf(", %d", number.modules_system[i]);
    }

    printf(")\n");
}

crt sum_crt_numbers(crt number1, crt number2) {
    crt result;

    result.problem_length = number1.problem_length;
    result.modules_system = number1.modules_system;
    result.modules_prod = number1.modules_prod;
    
    int *residuals = (int*) malloc(result.problem_length * sizeof(int));
    int i;
    for (i = 0; i < result.problem_length; i++) {
        residuals[i] = (number1.residuals[i] + number2.residuals[i]) % result.modules_system[i];
    }
    result.residuals = residuals;

    return result;
}


crt prod_crt_numbers(crt number1, crt number2) {
    crt result;
    
    result.problem_length = number1.problem_length;
    result.modules_system = number1.modules_system;
    result.modules_prod = number1.modules_prod;
    
    int *residuals = (int*) malloc(result.problem_length * sizeof(int));
    int i;
    for (i = 0; i < result.problem_length; i++) {
        residuals[i] = (number1.residuals[i] * number2.residuals[i]) % result.modules_system[i];
    }
    result.residuals = residuals;

    return result;
}

/*
Данная функция ищет обратное значение для числа по модулю 
в мультипликативной форме по правилу модулярной арифметики

@param int residual: остаток от деления числа по модулю
@param in module: модуль деления с остатком

@return int inverse: найденное обратное число по модулю

в случае если число нашлось, выдает его,
если нет то выдает ноль
*/
int inverse_residual(int residual, 
                     int module) {
    int inverse = 1;
    for (inverse = 1; inverse < module; inverse++) {
        if ((residual * inverse) % module == 1) 
            return inverse;
        if ((residual * inverse) % module == 0) 
            return 0;
    }
};


/*Данная функция решает задачу в постановке числа в виде СОК*/
int solve_china_resid(crt problem) {
    int result = 0;
    int *modules_primes = (int*) malloc(problem.problem_length * sizeof(int));
    int *inverse_array = (int*) malloc(problem.problem_length * sizeof(int));

    int i;
    for (i = 0; i < problem.problem_length; i++) {

        // N[i] -> n[1] * n[2] * ... * n[i-1] * n[i+1] * ... * n[k]
        modules_primes[i] = problem.modules_prod / problem.modules_system[i];
        
        // x[i] -> N[i] * x[i] = 1 {n[i]}
        inverse_array[i] = inverse_residual(modules_primes[i] % problem.modules_system[i], 
                                            problem.modules_system[i]);
        
        // += (x[i] * N[i] * r[i]) {n[i]}
        result += inverse_array[i] * modules_primes[i] * problem.residuals[i] % problem.modules_prod;
    
    }

    // X = X(l){N}; where X(l) = X + N * l, l = 0, 1, 2, ...
    return result % problem.modules_prod;
}


// Пример работы поиска обратных чисел по модулю в мультипликативной форме
int example_prog1() {
    int array_length = 5;
    int resid[5] = {5, 7, 9, 11, 13};
    int module[5] = {7, 9, 11, 13, 17};
    
    int arr_index; 
    for (arr_index = 0; arr_index < array_length; arr_index++) {
        printf("inverse of %d{%d} is %d{%d}\n", 
               resid[arr_index], 
               module[arr_index], 
               inverse_residual(resid[arr_index], module[arr_index]), 
               module[arr_index]);
    }
    return 0;
}


// Пример работы объявления структуры числа в СОК 
// постановка задачи КТО 
// решение задачи КТО
int example_prog2() {
    int array_length = 5;
    int resid[5] = {5, 7, 9, 11, 13};
    int module[5] = {7, 9, 11, 13, 17};
    
    // Пример работы объявления структуры числа в СОК
    crt sample_problem = init_crt(array_length, resid, module, 0);

    // постановка задачи КТО 
    print_crt_problem(sample_problem);

    // печать структуры числа в СОК
    print_crt_number(sample_problem);
    
    // решение задачи КТО
    printf("%d\n", solve_china_resid(sample_problem));
    return 0;
}


/*Пример полного цикла арифметики в СОК*/
int example_prog3() {
    int system_len = 4;      // Размерность системы модулей равна l = 4
    int modules_system[4] = {12, 7, 17, 5};  // Модули системы n = (n_1, n_2, ..., n_l)
    int module_base = prod_array_i(modules_system, system_len); // Произведение модулей N
    int A = 131, B = 15, C = 46;

    printf("Modules system is: n = (%d", modules_system[0]);
    int i;
    for (i = 1; i < system_len; i++) {
        printf(", %d", modules_system[i]);
    }
    
    printf(")\nNumbers is\n A = %d \n B = %d \n C = %d \n", A, B, C);
    printf("\nIn residuals of modular arifmetics numbers A, B, C is:\n");
    crt A_crt, B_crt, C_crt;

    A_crt = init_crt_from_int(system_len, A, modules_system, module_base);
    print_crt_number(A_crt);
    
    B_crt = init_crt_from_int(system_len, B, modules_system, module_base);
    print_crt_number(B_crt);
    
    C_crt = init_crt_from_int(system_len, C, modules_system, module_base);
    print_crt_number(C_crt);
    

    crt A_plus_B_crt = sum_crt_numbers(A_crt, B_crt);
    if (solve_china_resid(A_plus_B_crt) == A + B) {
        printf("\nA + B = Acrt + B_crt\n");
    }
    
    crt A_plus_C_crt = sum_crt_numbers(A_crt, C_crt);
    if (solve_china_resid(A_plus_C_crt) == A + C) {
        printf("A + C = Acrt + C_crt\n");
    }
    
    crt C_plus_B_crt = sum_crt_numbers(C_crt, B_crt);
    if (solve_china_resid(C_plus_B_crt) == B + C) {
        printf("B + C = Ccrt + B_crt\n");
    }

    crt A_time_B_crt = prod_crt_numbers(A_crt, B_crt);
    if (solve_china_resid(A_time_B_crt) == A * B) {
        printf("A * B = A_crt * B_crt\n");
    }
    
    crt A_time_C_crt = prod_crt_numbers(A_crt, C_crt);
    if (solve_china_resid(A_time_C_crt) == A * C) {
        printf("A * C = A_crt * C_crt\n");
    }

    crt C_time_B_crt = prod_crt_numbers(C_crt, B_crt); 
    if (solve_china_resid(C_time_B_crt) == C * B) {
        printf("B * C = B_crt * C_crt\n");
    }

    return 0;
}


int make_excersices_latex() {
    int system_length = 3;
    int examples = 30;
    int modular_systems[30][3] = {{4, 7, 29}, 
                                  {3, 11, 13},
                                  {5, 11, 17},
                                  {11, 19, 23},
                                  {5, 7, 19},
                                  {11, 14, 17},
                                  {3, 7, 29},
                                  {11, 12, 13},
                                  {5, 7, 29}, 
                                  {4, 10, 11},
                                  {5, 11, 17},
                                  {10, 19, 23},
                                  {5, 7, 24},
                                  {11, 14, 19},
                                  {3, 11, 29},
                                  {11, 12, 23},
                                  {5, 12, 29}, 
                                  {4, 11, 17},
                                  {5, 12, 17},
                                  {11, 19, 24},
                                  {5, 7, 11},
                                  {11, 14, 24},
                                  {3, 7, 20},
                                  {11, 12, 13},
                                  {5, 11, 29}, 
                                  {4, 11, 17},
                                  {3, 10, 17},
                                  {11, 13, 23},
                                  {5, 13, 19},
                                  {10, 11, 17}};
    int X = 97;
    int i;

    crt problem;
    for (i = 0; i < 30; i++) {
        problem = init_crt_from_int(3, X, modular_systems[i], prod_array_i(modular_systems[i], 3));
        print_crt_problem_latex(problem);
    }
    return 0;
}


int main() {
    //example_prog1();
    //example_prog2();
    //example_prog3();
    make_excersices_latex();
    return 0;
}