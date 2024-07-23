# Contains useful functions for parsing strings containing Boolean instructions, appying these
# conditions to a Pandas DataFrame, and returning a Pandas series of Booleans inidicating whether
# the conditions are met for each row

import operator
import re

def str_2_logical(in_string):
    # Convert string containing logical (or arithmetic) expression to evaluate the expression itself
    mapping = {"<"  : operator.lt,
               "<=" : operator.le,
               ">"  : operator.gt,
               ">=" : operator.ge,
               "&"  : operator.and_,
               "|"  : operator.or_,
               "+"  : operator.add,
               "-"  : operator.sub,
               "*"  : operator.mul,
               "/"  : operator.truediv,
               "~"  : operator.not_,
               "!=" : operator.ne}
    
    return(mapping[in_string])

def partition_by_op(in_string, ops):
    # Partition string according to specified operators
    out_list = [in_string]
    out_ops = []
    for op in ops:
        for i, string in enumerate(out_list):
            if op in string:
                # I've already added spaces around the minus sign to distinguish between minus and subtract
                if not op == " - ": 
                    # Pad with spaces if not already
                    out_list[i] = string.split(op.center(len(op)+2))
                else:
                    out_list[i] = string.split(op)
                out_list = [item for item in out_list]
                # Flatten the list to avoid nesting
                flt_list = []
                for item in out_list:
                    if isinstance(item, list):
                        flt_list.extend(item)
                    else:
                        flt_list.append(item)
                out_list = flt_list

    if len(out_list) == 1:
        out_list = out_list[0]
    
    # Use regular expressions to return a list of the operators in the order in which they appear
    for char in ["|","+","-","*"]:
        ops = [op.replace(char,"\\" + char) for op in ops]

    pattern = r"|".join(ops)
    out_ops = re.findall(pattern, in_string)

    return(out_list, out_ops)

def get_log_index(data, in_string):
    # Get logical index for applying a set of instructions, contained in a string, to a dataframe
    # Avoid using "eval" as this is bad practice for security reasons
    # Partition by and and or operations
    indiv_conds, linking_ops = partition_by_op(in_string, ["&","|"])
    # how deep is indiv_ops? If only one convert to list for looping
    if len(indiv_conds[0]) == 1:
        indiv_conds = [indiv_conds]

    indiv_inds = []
    for cond in indiv_conds:
        # Now partition string by inequalities
        cond, indiv_op = partition_by_op(cond, ["<=",">=",">","<"]) # Check less than or equal to first otherwise they get wrongly identified
        # If a simple condition the second term should be a float. If this doesn't work further partitioning is required
        try:
            cond[1] = float(cond[1])
            # Append logical index to list
            # If negating the values in the column remove the "-" from the index, use the rest of the string and negate the data
            if cond[0][0] == "-":
                indiv_inds.append(str_2_logical(indiv_op[0])(-data[cond[0][1:]],cond[1]))
            else:
                indiv_inds.append(str_2_logical(indiv_op[0])(data[cond[0]],cond[1]))
        except:
            # Further partition second term for arithmetic operations
            # Add spaces to " - " to distinguish subtract from minus
            complex_cond, complex_cond_ops = partition_by_op(cond[1], ["+"," - ","*","/"])
            # Extract the data from thte dataframe and apply arithmetic operations,
            # then compute the logical index
            # If negating the values in the column remove the "-" from the index, use the rest of the string and negate the data
            if complex_cond[0][0] == '-':
                data_cond = -data[complex_cond[0][1:]]
            else:
                data_cond = data[complex_cond[0]]
            numeric_vals = [float(item) for item in complex_cond[1:]]
            for op, val in zip(complex_cond_ops, numeric_vals):
                data_cond = str_2_logical(op)(data_cond,val)

            # If negating the values in the column remove the "-" from the index, use the rest of the string and negate the data
            if cond[0][0] == '-':
                indiv_inds.append(str_2_logical(indiv_op[0])(-data[cond[0][1:]],data_cond))
            else:
                indiv_inds.append(str_2_logical(indiv_op[0])(data[cond[0]],data_cond))
    
    # Get function for linking and/or term
    linking_ops = [str_2_logical(op) for op in linking_ops]
    # Link together logical indices with and/or as required
    out_ind = indiv_inds[0]
    if len(linking_ops) > 0:
        for op, link in zip(indiv_inds[1:],linking_ops):
            out_ind = link(out_ind, op)

    return(out_ind)