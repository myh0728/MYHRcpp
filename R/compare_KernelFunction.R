K2_Epanechnikov_ifelse <- function(u)
{
  ifelse(abs(u) <= 1, 3 / 4 * (1 - u ^ 2), 0)
}
