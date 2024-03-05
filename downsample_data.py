def downsample_data(data, sort = False, sort_col = ["Element"], method = "uniform", rate = 2):
    # Downsample the data to reduce sample size
    n_data = data.shape[0]
    # Sort into order of ascending element number such that points in close proximity are grouped
    if sort:
        data.sort_values(sort_col,inplace=True)
    if method == "uniform":
        # Take a uniform subset of the samples by selecting at regular intervals
        data = data.iloc[range(0,n_data,rate)]
    elif method == "random":
        # randomly select the desired subset of the data
        data = data.sample(frac = 1.0/rate)
    return(data)