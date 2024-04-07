states <- c("Rouen", "Bucharest", "Samos", "Aigio", "Marseille")
emc <- create_sequence(states, probs = c(0.3, 0.1, 0.1, 0.3, 0.2))
obj_model <- fit_dsmm(
    sequence = emc,
    states = states,
    degree = 3,
    f_is_drifting = FALSE,
    p_is_drifting = TRUE
)

test_that("get_kernel()", {

    kern <- get_kernel(obj_model, klim = 10)
    expect_type(kern, "double")
    expect_identical(
        # 3rd & 4th dimension will vary across runs
        dim(kern)[1:2],
        c(length(states), length(states))
    )
    # Probabilities sum to 1
    expect_equal(
        rowSums(kern, 1),
        setNames(rep_len(dim(kern)[4], length(states)), states)
    )

})

test_that("simulate()", {

    s <- simulate(obj_model, max_seq_length = 10)

    expect_lte(length(s), 10)
    expect_in(s, states)

})
