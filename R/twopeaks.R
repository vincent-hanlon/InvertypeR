#' Finds two peaks in a vector
#'
#' Looks for the two biggest peaks (local maxima) in a vector of non-negative integers.
#'
#' This is used to determine the start and end coordinates of an inversion from a list of deltaW values, which are associated with reads. It uses run length encoding
#' to find the two highest local maxima. It doesn't select two adjacent peaks: instead, it introduces a spacer equal to 1/12th the length of the vector. It returns the
#' coordinates of the right side of the leftmost peak and the left side of the rightmost peak.
#'
#' @param x A vector of non-negative integers.
#' @return A vector of length two.
#' @examples
#'
#' twopeaks(c(1, 2, 4, 2, 1, 3, 3, 3, 5, 5, 5, 6, 6, 2, 2, 2, 5, 4, 5, 4, 7))
twopeaks <- function(x) {

    t <- rle(x)

    if (length(t$length) > 1) {
        # retain only local maxima in t
        t$values[-(which(diff(sign(diff(t$values))) == -2) + 1)] <- 0

        # Finding the highest peak
        m1 <- which.max(t$values)

        # Setting a rough spacer: peaks shouldn't be too close together!
        space <- as.integer(length(t$values) / 12)

        # Set highest peak and nearby values to zero
        t$values[max(m1 - space, 1):min(m1 + space, length(t$values))] <- 0

        # Finding the second-highest peak, and figuring out which one comes first
        m2 <- which.max(t$values)
        m <- sort(c(m1, m2))

        # Returning the inner coordinates of the peaks


        coords <- c(sum(t$lengths[1:m[1]]), sum(t$lengths[1:(m[2] - 1)]) + 1)
    } else if (length(t$length) == 1) {
        coords <- (c(1, length(x)))
    }

    return(coords)
}
