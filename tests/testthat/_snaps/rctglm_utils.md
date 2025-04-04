# `oos_fitted.values_counterfactual` snapshot test

    Code
      oos
    Output
         psi0 psi1 rowname
      1     1    1       1
      2     2    2       2
      3     3    3       3
      4     4    4       4
      5     5    5       5
      6     6    6       6
      7     7    7       7
      8     8    8       8
      9     9    9       9
      10   10   10      10

# `extract_train_test` returns list of train and test data

    Code
      list_of_train_test
    Output
      $train
        Y X A
      3 3 3 0
      4 4 4 0
      5 5 5 0
      6 6 6 1
      7 7 7 1
      
      $test
          Y  X A
      1   1  1 0
      2   2  2 0
      8   8  8 1
      9   9  9 1
      10 10 10 1
      

