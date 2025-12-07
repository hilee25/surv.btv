#' mantel_byar_cp(): Mantel-Byar 검정 함수
#' 각 사건 시점에서 위험집단을 재구성하여 Mantel-Byar 검정을 수행
#' - 입력: counting-process 데이터(id, tstart, tstop, event, response), ref_val(기준 상태 값)
#' - 처리: 사건 발생 시점들의 위험집단에서 계산
#' - 출력: list(chisq, p, D_ref, E, V, contrib(시점별 표), ref_val)
#'
#' @param cp data.frame. columns(id, tstart, tstop, event, response).
#' @param ref_val 0 또는 1. 통상 0을 참조군으로 둠.
#'
#' @return list(chisq, p, D_ref, E, V, contrib, ref_val)
#' @export
mantel_byar_cp <- function(cp, ref_val = 0L) {
  req <- c("tstart","tstop","event", "response")
  stopifnot(all(req %in% names(cp)))
  if (!is.numeric(cp$tstart) || !is.numeric(cp$tstop)) {
    stop("start/stop must be numeric times (e.g., days since baseline).")
  }
  if (any(cp$tstop <= cp$tstart, na.rm = TRUE)) {
    stop("Found non-positive intervals: require start < stop for all rows.")
  }

  z <- cp[["response"]]
  if (is.factor(z)) z <- as.character(z)
  z <- as.integer(z)
  badz <- which(!is.na(z) & !(z %in% c(0L,1L)))
  if (length(badz)) stop("response must be binary 0/1; offending rows: ", paste(head(badz,5), collapse=", "))
  cp[["response"]] <- z

  evt_times <- sort(unique(cp$tstop[cp$event == 1]))
  if (length(evt_times) == 0) {
    return(list(chisq = 0, p = 1, D_ref = 0, E = 0, V = 0,
                contrib = data.frame(time = numeric(0), d = integer(0),
                                     N = integer(0), Nref = integer(0),
                                     dref = integer(0), OE = numeric(0), V = numeric(0))))
  }

  D_ref <- 0
  E_sum <- 0
  V_sum <- 0

  contrib <- vector("list", length(evt_times))

  for (i in seq_along(evt_times)) {
    tt <- evt_times[i]

    at_risk <- cp[cp$tstart < tt & cp$tstop >= tt, , drop = FALSE]
    if (nrow(at_risk) == 0) {
      contrib[[i]] <- data.frame(time = tt, d = 0L, N = 0L, Nref = 0L,
                                 dref = 0L, OE = 0, V = 0)
      next
    }

    N    <- nrow(at_risk)
    Nref <- sum(at_risk[["response"]] == ref_val, na.rm = TRUE)

    at_evt <- cp[cp$tstop == tt & cp$event == 1L, , drop = FALSE]
    d      <- nrow(at_evt)
    if (d == 0L) {
      contrib[[i]] <- data.frame(time = tt, d = 0L, N = N, Nref = Nref,
                                 dref = 0L, OE = 0, V = 0)
      next
    }
    dref   <- sum(at_evt[["response"]] == ref_val, na.rm = TRUE)

    if (N <= 1L || Nref == 0L || Nref == N) {
      contrib[[i]] <- data.frame(time = tt, d = d, N = N, Nref = Nref,
                                 dref = dref, OE = 0, V = 0)
      next
    }

    E  <- d * (Nref / N)
    V  <- d * (Nref * (N - Nref) / (N^2)) * ((N - d) / (N - 1))

    OE <- (dref - E)

    D_ref <- D_ref + dref
    E_sum <- E_sum + E
    V_sum <- V_sum + V

    contrib[[i]] <- data.frame(time = tt, d = d, N = N, Nref = Nref,
                               dref = dref, OE = OE, V = V)
  }

  contrib <- do.call(rbind, contrib)

  if (!is.finite(V_sum) || V_sum <= 0) {
    chisq <- 0
    pval  <- 1
  } else {
    chisq <- ( (D_ref - E_sum)^2 ) / V_sum
    pval  <- pchisq(chisq, df = 1, lower.tail = FALSE)
  }

  list(chisq = chisq, p = pval, D_ref = D_ref, E = E_sum, V = V_sum,
       contrib = contrib, ref_val = ref_val)
}
