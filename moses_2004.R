get_delta <- function(k, l, rho) {
  numerator = (k/l)^(1+rho)
  return(numerator/(1 + numerator))
}

get_gamma <- function(gdp, delta, k, l, rho) {
  denominator = (delta * k^(-rho) + (1-delta)*l^(-rho))^(1/-rho)
  return(gdp/denominator)
}

get_mpl <- function(gdp, delta, gamma, l, rho) {
  return((1-delta) * gamma^(-rho) * gdp^(1+rho) * l^(-rho-1))
}

get_gdp <- function(delta, gamma, k, l, rho) {
  return(gamma * (delta * k^(-rho) + (1-delta) * l^(-rho))^(1/-rho))
}

get_labor <- function(mpl, gdp, delta, gamma, rho) {
  return((mpl / ((1-delta) * gamma^(-rho) * gdp^(1+rho)))^(1/(-rho-1)))
}

get_rho <- function(sigma) {
  return(1/sigma - 1)
}

get_k <- function(gdp, l) {
  # labour share, based on UN National Accounts Statistics (1986, 1997 & 2000)
  ls=.687;
  w = ls * gdp / l
  return((1-ls)*gdp/w) # gdp is not even necessary
}

compute <- function(workforce, gdp) {
  capital <- get_k(gdp, workforce)

  partial_delta <- function(rho) {
    return(get_delta(capital, workforce, rho))
  }
  deltas <- data.frame(sapply(rhos, partial_delta))
  colnames(deltas) <- rhos # index deltas columns by rho for lookup later

  partial_gamma <- function(rho) {
    return(get_gamma(gdp, deltas[as.character(rho)],
                     capital, workforce, rho))
  }
  gammas <- data.frame(sapply(rhos, partial_gamma))
  colnames(gammas) <- rhos

  try_solve_for_l <- function(i, rho, mpl) {
    deltai <- deltas[i, as.character(rho)]
    gammai <- gammas[i, as.character(rho)]
    ki <- capital[i, 1]
    fn <- function(li) {
      gdpi <- get_gdp(deltai, gammai, ki, li, rho)
      return(get_labor(mpl, gdpi, deltai, gammai, rho) - li)
    }
    labor <- NA
    try(labor <- uniroot(fn, c(1e-20, total_pop), tol=1e-9)$root, # accurate to 1 person, Pop in billions
      silent=T)
    if(is.na(labor)) {
      return(c(NA, NA))
    } else {
      return(c(labor, get_gdp(deltai, gammai, ki, labor, rho)))
    }
  }

  INITIAL_MIN = 100
  INITIAL_MAX = 50000
  find_mpl_range <- function(rho) {
    min = INITIAL_MIN
    max = INITIAL_MAX
    for (i in 1:3) {
      for (mpl in seq(min, max, 100)) {
        result <- try_solve_for_l(i, rho, mpl)
        if (!anyNA(result)) {
          min = mpl
          break
        }
      }
      # switch to smaller step size
      for (mpl in seq(min-100, min, 1)) {
        result <- try_solve_for_l(i, rho, mpl)
        if (!anyNA(result)) {
          min = mpl
          break
        }
      }

      for (mpl in seq(max, min, -100)) {
        result <- try_solve_for_l(i, rho, mpl)
        if (!anyNA(result)) {
          max = mpl
          break
        }
      }
      # switch to smaller step size
      for (mpl in seq(max+100, max, -1)) {
        result <- try_solve_for_l(i, rho, mpl)
        if (!anyNA(result)) {
          max = mpl
          break
        }
      }
      if (min != INITIAL_MIN && max != INITIAL_MAX) {
        break
      }
    }

    return(c(min, max))
  }

  # simulate immigration liberalization
  # search through a range of MPLs
  total_pop = sum(workforce)
  for (i in 1:length(rhos)) {
    min_pop_diff = .Machine$integer.max
    rho = rhos[i]
    solved <- data.frame()
    solved_mpl <- NaN
    mpl_range <- find_mpl_range(rho)
    # print(mpl_range)
    for (mpl in seq(mpl_range[1], mpl_range[2], 10)) {
      # solve for l1, l2 and l3 with the given mpl
      solve_for_l <- function(j) {
        return(try_solve_for_l(j, rho, mpl))
      }

      new_labor_gdp <- data.frame(t(sapply(1:3, solve_for_l)))
      colnames(new_labor_gdp) <- c('labor', 'gdp')

      total_l = sum(new_labor_gdp['labor'], na.rm = T)
      pop_diff = abs(total_l - total_pop)
      if (pop_diff < min_pop_diff) {
        min_pop_diff = pop_diff
        solved = new_labor_gdp
        solved_mpl = mpl
      }
    }
    print(solved)
    print(paste('sigma:', elasticities[i], 'rho:', rho, ' mpl:', solved_mpl))
    print(paste('abs pop diff (billion):', min_pop_diff))
    print(paste('gdp diff (billion $):', sum(solved['gdp'], na.rm = T) - sum(gdp)))
    print('-------------------------------------------')
  }
}

pw_adjust <- function(workforce) {
  return(c(workforce[1]*0.48, workforce[2]*0.41, workforce[3]*0.41))
}

eu35_adjust <- function(workforce) {
  return(c(workforce[1], workforce[2]/3, workforce[3]/5))
}

elasticities = c(0.5, 0.75, 0.999, 1.25, 1.5)
rhos = get_rho(elasticities)

pop_gdp <- read.csv("borderstest2.csv")
gdp <-pop_gdp['GDP']

# no adjustment
workforce <- pop_gdp['Pop']
compute(workforce, gdp)

# PW adjustment
compute(pw_adjust(workforce[, 1]), gdp)

# EU3 & EU5 adjustment
compute(eu35_adjust(workforce[, 1]), gdp)

# both adjustments
compute(eu35_adjust(pw_adjust(workforce[, 1])), gdp)


