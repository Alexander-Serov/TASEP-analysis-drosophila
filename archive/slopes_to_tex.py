

from constants import (slope_abortive_theory, slope_abortive_theory_error,
                       tex_slopes_data_file)


def slopes_to_tex(slopes):
    """
    Output slope values into a tex file to use for the article
    """

    def add_string(tex_file, nc, column_label, theor_value, theor_str, head_str=""):
        out_str = head_str
        out_str += "\t&\t%i" % (nc)
        for gene_id in range(3):
            for construct_id in range(3):

                # Calculate current slope
                cur_data = slopes[(slopes.gene_id == gene_id) & (
                    slopes.construct_id == construct_id)][column_label].dropna()

                mean = cur_data.mean()
                std = np.std(cur_data, ddof=1)

                # Prepare output string
                out_str += "\t&\t$%.0f\\pm%.0f$" % (mean, std)

        # Write to file
        tex_file.write(out_str)

    with open(tex_slopes_data_file, "w") as tex_file:

        # s, nc 13
        nc = 13
        theor_str = "$%.0f\pm%.0f$" % (slope_abortive_theory, slope_abortive_theory_error)
        head_str = "\\multirow{2}{\\multirowWidth}{$s$}"
        theor_value = slope_abortive_theory
        add_string(tex_file, nc, 'slope_nc13', theor_value, theor_str, head_str)
        tex_file.write("\n\\\\\n")

        # s, nc 14
        nc = 14
        add_string(tex_file, nc, 'slope_nc14', theor_value, theor_str)
        tex_file.write("\n\\vspace{2mm}\n\\\\\n")

        # N, nc 13
        nc = 13
        head_str = "\\multirow{2}{\\multirowWidth}{$\\NSS$}"
        theor_str = ""
        theor_value = 1
        add_string(tex_file, nc, 'max_nc13', theor_value, theor_str, head_str)
        tex_file.write("\n\\\\\n")

        # N, nc 14
        nc = 14
        add_string(tex_file, nc, 'max_nc14', theor_value, theor_str)
        tex_file.write("\n\\\\\n")
