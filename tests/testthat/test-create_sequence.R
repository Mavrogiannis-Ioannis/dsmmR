test_that("create_sequence() gives expected names and frequencies", {

    # This is equal to having the argument `probs = c(1/4, 1/4, 1/4, 1/4)`.
    set.seed(1)
    rand_dna_seq <- create_sequence(states = "DNA", len = 1e7, seed = 1)
    expect_equal(
        c(table(rand_dna_seq)) / 1e7,
        setNames(rep_len(0.25, 4), c("a", "c", "g", "t")),
        tolerance = 1e-2
    )

    set.seed(1)
    probs <- c(0.6, 0.3, 0.05, 0.025, 0.025)
    random_letters <- sample(letters, size = 5, replace = FALSE)
    rand_dna_seq2 <- create_sequence(
        states = random_letters,
        probs = probs,
        len = 1e7,
        seed = 1
    )

    expected <- setNames(probs, unique(random_letters))
    expected <- expected[order(names(expected))]

    expect_equal(
        c(table(rand_dna_seq2)) / 1e7,
        expected,
        tolerance = 1e-2
    )

})
