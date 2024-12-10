import numpy as np
from numpy.random import randint
from sklearn.preprocessing import OneHotEncoder



def get_mask(view_num, alldata_len, missing_rate):
    """Randomly generate incomplete data information, simulate partial view data with complete view data
    :param view_num:view number
    :param alldata_len:number of samples
    :param missing_rate:Defined in section 4.1 of the paper
    :return: mask
    """
    full_matrix = np.ones((int(alldata_len * (1 - missing_rate)), view_num))

    alldata_len = alldata_len - int(alldata_len * (1 - missing_rate))
    missing_rate = 0.5
    if alldata_len != 0:
        one_rate = 1.0 - missing_rate
        if one_rate <= (1 / view_num):
            enc = OneHotEncoder()  # n_values=view_num
            view_preserve = enc.fit_transform(randint(0, view_num, size=(alldata_len, 1))).toarray()
            full_matrix = np.concatenate([view_preserve, full_matrix], axis=0)
            choice = np.random.choice(full_matrix.shape[0], size=full_matrix.shape[0], replace=False)
            matrix = full_matrix[choice]
            return matrix
        error = 1
        if one_rate == 1:
            matrix = randint(1, 2, size=(alldata_len, view_num))
            full_matrix = np.concatenate([matrix, full_matrix], axis=0)
            choice = np.random.choice(full_matrix.shape[0], size=full_matrix.shape[0], replace=False)
            matrix = full_matrix[choice]
            return matrix
        while error >= 0.005:
            full_range = np.arange(view_num)
            shuffled_range = np.random.permutation(full_range)
            # Repeat the shuffled array until it reaches the required length
            randint_array = np.tile(shuffled_range, int(np.ceil(alldata_len / view_num)))
            # Trim the array to the exact length needed
            randint_array = randint_array[:alldata_len].reshape(-1, 1)

            enc = OneHotEncoder()
            view_preserve = enc.fit_transform(randint_array).toarray()

            # enc = OneHotEncoder()  # n_values=view_num
            # view_preserve = enc.fit_transform(randint(0, view_num, size=(alldata_len, 1))).toarray()
            # print(view_preserve.shape)
            # if view_preserve.shape[1]==3:
            #     view_preserve = enc.fit_transform(randint(0, view_num+1, size=(alldata_len, 1))).toarray()
            one_num = view_num * alldata_len * one_rate - alldata_len
            ratio = one_num / (view_num * alldata_len)
            matrix_iter = (randint(0, 100, size=(alldata_len, view_num)) < int(ratio * 100)).astype(np.int32)
            # print(matrix_iter.shape)
            a = np.sum(((matrix_iter + view_preserve) > 1).astype(np.int32))
            one_num_iter = one_num / (1 - a / one_num)
            ratio = one_num_iter / (view_num * alldata_len)
            matrix_iter = (randint(0, 100, size=(alldata_len, view_num)) < int(ratio * 100)).astype(np.int32)
            matrix = ((matrix_iter + view_preserve) > 0).astype(np.int32)
            ratio = np.sum(matrix) / (view_num * alldata_len)
            error = abs(one_rate - ratio)
        full_matrix = np.concatenate([matrix, full_matrix], axis=0)

    choice = np.random.choice(full_matrix.shape[0], size=full_matrix.shape[0], replace=False)
    matrix = full_matrix[choice]
    return matrix

