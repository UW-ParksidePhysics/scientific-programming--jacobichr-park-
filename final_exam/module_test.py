modules = ['final_exam/read_two_columns_text', 'final_exam/calculate_bivariate_statistics', 'final_exam/calculate_quadratic_fit', 'final_exam/fit_curve_array', 'final_exam/plot_data_with_fit', 'final_exam/calculate_lowest_eigenvectors', 'final_exam/annotate_plot']
for module in modules:
    with open(module+'.py') as module_file:
        exec(module_file.read())