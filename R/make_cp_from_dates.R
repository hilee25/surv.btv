#' make_cp_from_dates(): Counting process 자료 변환 함수
#' 날짜 기반 원자료 → counting-process + time-dependent binary covariate 생성
#' - 입력: landmark_days, accept, response, follow-up end, status
#' - 출력: id, tstart, tstop, event, response 가 있는 data.frame
#'
#' @param data            분석할 data.frame
#' @param landmark_days   numeric 벡터. baseline(accept)로부터 며칠 후를 landmark로 지정 (예: c(60, 120)).
#' @param accept_col      character. baseline 날짜 변수명 (예: "accept.dt").
#' @param tx_col          character. 반응 날짜 변수명 (예: "tx.date").
#' @param fu_col          character. 추적 종료(사망·중도절단 발생 시점) 날짜 변수명 (예: "fu.date").
#' @param status_col      character. 사건 지표 변수명 (예: "fustat", 1=사망, 0=생존).
#' @param eps             0-길이 구간을 피하기 위해 아주 작은 시간을 더해 주는 보정값.
#' @param fix_zero_followup futime==0(accept=fu 같은 날) 처리 방식(bump: 아주 살짝 늘림, drop: 해당 사례 제외, error: 중단)
#'
#' @return data.frame. columns(id, tstart, tstop, event, response).
#' @export
make_cp_from_dates <- function(data,
                               accept_col,
                               tx_col,
                               fu_col,
                               status_col,
                               eps = 1e-8,
                               fix_zero_followup = c("bump","drop","error")) {

  fix_zero_followup <- match.arg(fix_zero_followup)

  stopifnot(is.data.frame(data))

  req <- c(accept_col, tx_col, fu_col, status_col)
  if (!all(req %in% names(data))) {
    stop("필수 열 누락: ", paste(setdiff(req, names(data)), collapse = ", "))
  }

  if (!("id" %in% names(data))) data$id <- seq_len(nrow(data))

  acc <- as.Date(data[[accept_col]])
  tx  <- as.Date(data[[tx_col]])
  fu  <- as.Date(data[[fu_col]])

  st <- data[[status_col]]
  if (is.factor(st)) st <- as.integer(as.character(st))
  st <- as.integer(st)

  futime <- as.numeric(fu - acc) # 추적기간
  txtime <- as.numeric(tx - acc)  # 반응까지의 시간 (NA 허용: 반응 없음)


  bad <- which(futime < 0)
  if (length(bad)) {
    stop("fu_col < accept_col 인 레코드가 있습니다. id: ",
         paste(head(data$id[bad], 5), collapse = ", "),
         if (length(bad) > 5) ", ...", "")
  }

  zero_fu <- which(futime == 0 | is.na(futime))
  if (length(zero_fu)) {
    if (fix_zero_followup == "bump") {
      futime[zero_fu] <- futime[zero_fu] + eps
    } else if (fix_zero_followup == "drop") {
      keep <- setdiff(seq_len(nrow(data)), zero_fu)
      data <- data[keep, , drop = FALSE]
      acc  <- acc[keep]; tx <- tx[keep]; fu <- fu[keep]
      st   <- st[keep];  futime <- futime[keep]; txtime <- txtime[keep]
    } else {
      stop(sprintf("추적 0일 사례 %d건 error 발생", length(zero_fu)))
    }
  }

  txtime[!is.na(txtime) & txtime < 0] <- NA_real_
  txtime[!is.na(txtime) & txtime == 0] <- eps

  idx_ge <- which(!is.na(txtime) & !is.na(futime) & txtime >= futime)
  if (length(idx_ge)) {
    txtime[idx_ge] <- NA_real_
  }

  base <- data.frame(id = data$id, tstart = 0, tstop = futime)

  tm <- tmerge(data1 = base, data2 = base, id = id,
               event = event(tstop, st == 1L))

  tm <- tmerge(tm, base, id = id,
               response = tdc(txtime))

  tiny <- 1e-12
  tm <- tm[(tm$tstop - tm$tstart) > tiny, , drop = FALSE]

  if (any(tm$tstop <= tm$tstart)) {
    stop("여전히 0-길이 구간이 남아 있습니다. eps를 키우세요.")
  }

  tm[, c("id","tstart","tstop","event","response")]
}
