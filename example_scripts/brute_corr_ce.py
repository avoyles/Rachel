
# The methods used in this script are not simple user accessors, so the script
# may be difficult to understand or modify.  For help with scripting, please
# make a post on the Gosia forums:
# http://www-user.pas.rochester.edu/~gosia/phpBB3


print "This version calculates new correction factors at each step for better accuracy.\nIt takes much longer than brute_ce.py."

if not yes_no_prompt("Have you done at least one chi-squared fit? [Y/n]: ",True):
    print "You must have completed at least one iteration of minimization (fitting) WITH THE CURRENT FIT PARAMETERS before running this script.  Quitting."

else:
    # The script can be run.

    # Save a copy of the original matrix.
    original_matrix = copy.deepcopy(investigated_nucleus.matrix_data)

    block_print_with_line_breaks("This script calculates the correlated error in one parameter at a time, that is, one master matrix element that you have already defined and all of its dependents as one parameter.\nTo calculate the full set of correlated errors treating each master and dependent as a *separate* fit parameter, use the \"Correlated errors\" function, option \"r\".\nThis script includes only correlations with other adjustable parameters, not with fixed matrix elements.\nThe purpose of this script is to provide a second check of the parametric option (\"p\") in the Gosia calculations, since this is not a standard Gosia operation.\nNote also that this script does a fast set of calculations that saves time by using the corrected yields near the minimum.  This will favor values nearest to the minimum and may underestimate the error.  You should compare this output to the built-in correlated error calculation.")

    block_print_with_line_breaks("Please choose a matrix element number from this list.  The correlated error will be calculated for this matrix element coupled to its dependents, if there are any.")

    # Show all master matrix elements, and prompt for the one to calculate the error for.
    me_number = 1
    all_matrix_keys = investigated_nucleus.matrix_data.keys()
    all_matrix_keys.sort()
    master_matrix_keys = []
    for one_matrix_key in all_matrix_keys:
        if investigated_nucleus.matrix_data[one_matrix_key].get_is_master():
            multipole_code, initial_band_name, initial_spin, final_band_name, final_spin = one_matrix_key
            multipole_text = REVERSE_MULTIPOLE[multipole_code]
            print str(me_number).rjust(3) + " " + investigated_nucleus.format_one_matrix_element(initial_band_name,initial_spin,final_band_name,final_spin,multipole_text)
            master_matrix_keys.append(one_matrix_key)
            me_number += 1

    # Get the number of degrees of freedom.  Subract the fixed one, because
    # this is how Gosia calculates the reduced chi-squared.
    degrees_of_freedom = len(master_matrix_keys) - 1

    chosen_me_number = int(raw_input("Matrix element number? "))

    chosen_matrix_key = master_matrix_keys[chosen_me_number - 1]
    chosen_matrix_object = investigated_nucleus.matrix_data[chosen_matrix_key]
    original_value = chosen_matrix_object.get_current_value()

    print "Current value: " + str(original_value)

    # Get the limits on this matrix element.
    lower_limit = float(raw_input("Enter the lower limit for the search: "))
    upper_limit = float(raw_input("Enter the upper limit for the search: "))
    number_of_sample_points = int(raw_input("Enter the number of sample points between these limits: "))
    print "If the correlated error cannot be found, broaden the limits."
    # Generate sample points for the master.
    master_sample_points = numpy.linspace(lower_limit,upper_limit,number_of_sample_points)

    # Find all dependents, and set them as fixed, so the fits won't change them.
    all_chosen_matrix_objects = []
    for this_matrix_object in investigated_nucleus.matrix_data.values():
        if this_matrix_object.get_master_matrix_element_key() == chosen_matrix_key:
            all_chosen_matrix_objects.append(this_matrix_object)
            this_matrix_object.set_is_fixed()
    # Set the master as fixed.
    chosen_matrix_object.set_is_fixed()

    # A list to record chi-squared:
    chi_squared_values = []

    # Step through values of the matrix element and do a fit at each point.
    previous_master_value = original_value
    iteration = 1
    for this_master_value in master_sample_points:
        # Change the chosen m.e. values in proportion to this change in the master.
        this_factor = this_master_value / previous_master_value 
        chosen_matrix_object.set_current_value(this_master_value)
        for one_dependent_object in all_chosen_matrix_objects:
            temporary_value = one_dependent_object.get_current_value()
            new_value = temporary_value * this_factor
            one_dependent_object.set_current_value(new_value)

        # Make new corrected yields.
        gosia_shell_output = the_gosia_shell.generate("Make corrected yields","Run gosia input")
        # Do the minimization for the remaining free parameters.
        gosia_shell_output = the_gosia_shell.generate("Fit","Run gosia input")

        # Record this chi-squared value.
        this_reduced_chi_squared = the_gosia_shell.return_final_chi_squared()
        this_chi_squared = this_reduced_chi_squared * float(degrees_of_freedom)
        chi_squared_values.append(this_chi_squared)

        # Print progress.
        print "\n\nIteration " + str(iteration) + ": " + " master value: " + str(this_master_value) + " reduced chi-squared: " + str(this_reduced_chi_squared) + "\n\n"
        sys.stdout.flush()
        time.sleep(1)

        # Increment things for the next iteration.
        previous_master_value = this_master_value
        iteration += 1


    print "Sampled *total* chi-squared (not reduced):"
    print "------------------------------------------"
    for i in range(len(chi_squared_values)):
        print master_sample_points[i], chi_squared_values[i]

    block_print_with_line_breaks("\nYou can examine this plot and make a polynomial fit to find the best value (minimum chi-squared) and the correlated errors where the curve crosses the chi-squared + 1 line, or another criterion.")

    print "The original matrix was restored (not updated from this calculation)."
    investigated_nucleus.matrix_data = copy.deepcopy(original_matrix)

    # Make the plot.
    horizontal = [[master_sample_points[0],master_sample_points[-1]],[min(chi_squared_values) + 1.,min(chi_squared_values) + 1.]]
    plot_dictionary = {"":[master_sample_points,chi_squared_values],"min + 1":horizontal}
    quick_plot_n_sets(plot_dictionary)





