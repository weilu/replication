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

try_solve_for_l <- function(i, rho, mpl) {
  deltai <- deltas[i, as.character(rho)]
  gammai <- gammas[i, as.character(rho)]
  ki <- pop_gdp[i, 'Capital']
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

find_mpl_range <- function(rho) {
  min = 100
  max = 50000
  for (i in 1:3) {
    for (mpl in seq(min, max, 100)) {
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
  }

  return(c(min, max))
}

pop_gdp <- read.csv("borderstest2.csv")
pop_gdp['Capital'] <- get_k(pop_gdp['GDP'], pop_gdp['Pop'])

# elasticities = c(.999)
elasticities = c(0.5, 0.75, 0.999, 1.25, 1.5)
rhos = get_rho(elasticities)

partial_delta <- function(rho) {
  return(get_delta(pop_gdp['Capital'], pop_gdp['Pop'], rho))
}
deltas <- data.frame(sapply(rhos, partial_delta))
colnames(deltas) <- rhos # index deltas columns by rho for lookup later

partial_gamma <- function(rho) {
  return(get_gamma(pop_gdp['GDP'], deltas[as.character(rho)],
                   pop_gdp['Capital'], pop_gdp['Pop'], rho))
}
gammas <- data.frame(sapply(rhos, partial_gamma))
colnames(gammas) <- rhos

# simulate immigration liberalization
# search through a range of MPLs
total_pop = sum(pop_gdp['Pop'])
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
  print(pop_gdp)
  print(solved)
  print(paste('sigma:', elasticities[i], 'rho:', rho, ' mpl:', solved_mpl))
  print(paste('abs pop diff (billion):', min_pop_diff))
  print(paste('gdp diff (billion $):', sum(solved['gdp'], na.rm = T) - sum(pop_gdp['GDP'])))
  print('-------------------------------------------')
}


