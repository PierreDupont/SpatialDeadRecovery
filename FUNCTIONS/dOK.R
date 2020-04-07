dOK <- nimbleFunction(
  run = function(x = double(0), sxy = double(1), lower = double(1), upper = double(1), habitat = double(2), log = double()) {
    returnType(double())
    if(sxy[1] < lower[1]) return(-Inf)   # x-coordinates
    if(sxy[1] > upper[1]) return(-Inf)   # x-coordinates
    if(sxy[2] < lower[2]) return(-Inf)   # y-coordinates
    if(sxy[2] > upper[2]) return(-Inf)   # y-coordinates
    ## habitat constraint:
    if(habitat[trunc(sxy[2])+1, trunc(sxy[1])+1] == 0) return(-Inf)
    return(0)
  }
)


rOK <- nimbleFunction(
  run = function(n = integer(), sxy = double(1), lower = double(1), upper = double(1), habitat = double(2)) {
    print('this should never happen')
    returnType(double(0))
    return(0)
  }
)

registerDistributions(list(
  dOK = list(
    BUGSdist = 'dOK(sxy, lower, upper, habitat)',
    types = c('value = double(0)', 'sxy = double(1)', 'lower = double(1)', 'upper = double(1)', 'habitat = double(2)'),
    mixedSizes = TRUE
  )
))