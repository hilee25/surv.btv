#' .survfit_at_L(): 시각화 도우미
#' 랜드마크 시점 L에서 KM 생존곡선 적합을 위한 데이터 구성 + survfit 객체 생성
#' - 입력: counting-process 데이터 (id, tstart, tstop, event, response), L
#' - 처리: (1) L 이전 사건 제거, (2) L 이후 추적 있는 id만 유지,
#'         (3) L을 덮는 구간의 response로 그룹(Z_L) 고정, (4) KM 적합
#' - 출력: list(sf = survfit 또는 NULL, postL = 분석 데이터, ztab =  (id 기준) 그룹 table)
#'
#' @export
.survfit_at_L <- function(cp_df, L) {
  stopifnot(all(c("id","tstart","tstop","event","response") %in% names(cp_df)))

  had_event_before_L <- cp_df %>%
    group_by(id) %>% summarize(ev_before = any(event == 1 & tstop <= L), .groups="drop")
  at_risk_ids <- cp_df %>%
    group_by(id) %>% summarize(has_followup = any(tstop > L), .groups="drop")

  keep_ids <- had_event_before_L %>%
    inner_join(at_risk_ids, by = "id") %>% filter(!ev_before & has_followup) %>% pull(id)
  df_risk <- cp_df %>% filter(id %in% keep_ids)

  zL <- df_risk %>%
    filter(tstart <= L, tstop > L) %>% select(id, response) %>% distinct() %>% rename(Z_L = response)

  postL <- df_risk %>%
    filter(tstop > L) %>%
    mutate(tstart = pmax(tstart, L)) %>%
    filter(tstop > tstart) %>%
    left_join(zL, by = "id") %>%
    mutate(Z_L = factor(Z_L, levels = c(0,1), labels = c("Z=0","Z=1")))

  gtbl <- table(distinct(postL, id, Z_L)$Z_L)

  if (length(gtbl) < 2 || any(gtbl == 0)) {
    return(list(sf = NULL, postL = postL, ztab = gtbl))
  }
  sf <- survfit(Surv(tstart, tstop, event) ~ Z_L, data = postL)

  list(sf = sf, postL = postL, ztab = gtbl)
}


#' plot_all_landmarks(): 다중 랜드마크 민감도 분석 및 시각화 함수
#' 여러 랜드마크 L들에 대해: (1) fit_for_L 결과/요약표 생성 + (2) KM 곡선/위험표 패널 그리드
#' - 입력: cp_df, landmarks, 시각화 옵션
#' - 처리: 각 L마다 fit_for_L + .survfit_at_L 수행
#' - 출력: list(plot = patchwork grid, table = tibble 요약표, per_L = 내부 객체 목록)
#'
#' @param cp_df            data.frame(id, tstart, tstop, event, response).
#' @param landmarks        numeric 또는 길이 1 이상의 벡터.
#' @param robust           cluster(id) 기반 강건 분산 사용 여부(기본 TRUE).
#' @param conf_int         KM 곡선의 신뢰구간 표시 여부(기본 TRUE).
#' @param tmax             모든 패널에서 공통 x축 상한(미지정 시 자동).
#' @param palette          라인 팔레트 (길이 2 권장).
#' @param theme_fn ggplot2 테마 함수. (기본 ggplot2::theme_bw).
#' @param base_size        메인 생존곡선 패널의 기본 글자 크기.
#' @param legend_position   범례 위치(기본 top).
#' @param ncol patchwork   그리드에서 한 행에 배치할 패널 수(기본 2).
#' @param title_prefix     패널 제목에 사용할 접두 문자열.
#' @param show_titles      각 랜드마크 패널에 제목 표시 여부(기본 TRUE).
#' @param show_risktable   각 패널 하단에 위험표 표시 여부(기본 TRUE).
#' @param risk_table_height 위험표가 패널 전체 높이에서 차지할 비율.
#' @param risk_table_base_size 위험표의 기본 글자 크기.
#' @param annotate_hrci    HR과 95% CI를 주석으로 표시할지 여부(기본 TRUE).
#' @param annotate_mb_p    Mantel–Byar p를 주석으로 표시할지 여부(기본 TRUE).
#' @param annotate_cox_p   Cox score p를 주석으로 표시할지 여부(기본 FALSE).
#' @param ref_label        참조군으로 간주할 Z 레이블 (예: "Z=0").
#'
#' @return list(plot = patchwork plot, table = tibble 요약표, per_L = 각 L의 객체 리스트)
#' @export
plot_all_landmarks <- function(cp_df,
                               landmarks,
                               robust = TRUE,
                               conf_int = TRUE,
                               tmax = NULL,
                               palette = c("#1f77b4", "#d62728"),
                               theme_fn = ggplot2::theme_bw,
                               base_size = 12,
                               legend_position = "top",
                               ncol = 2,
                               title_prefix = "Landmark L = ",
                               show_titles = TRUE,
                               show_risktable = TRUE,
                               risk_table_height = 0.35,
                               risk_table_base_size = 10,
                               annotate_hrci = TRUE,
                               annotate_mb_p = TRUE,
                               annotate_cox_p = FALSE,
                               ref_label = "Z=0"
) {
  if (length(palette) < 2) palette <- rep_len(palette, 2)
  L_vec <- sort(unique(as.numeric(landmarks)))
  if (any(!is.finite(L_vec))) stop("landmarks must be finite numeric.")

  fits <- purrr::map(L_vec, ~tryCatch(
    fit_for_L(cp_df, .x, robust = robust),
    error = function(e) list(error = e, L = .x)
  ))
  sfs  <- purrr::map(L_vec, ~.survfit_at_L(cp_df, .x))

  tab <- purrr::map2_dfr(fits, sfs, function(fit, sflist) {
    if (!is.null(fit$error)) {
      tibble(L = fit$L, n_total = NA_integer_, n_events = NA_integer_,
             n_Z0 = NA_integer_, n_Z1 = NA_integer_,
             HR = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
             p_mb = NA_real_, p_cox = NA_real_,
             note = paste0("fit error: ", conditionMessage(fit$error)))
    } else {
      z0 <- unname(if ("Z=0" %in% names(fit$n_group)) fit$n_group[["Z=0"]] else NA_integer_)
      z1 <- unname(if ("Z=1" %in% names(fit$n_group)) fit$n_group[["Z=1"]] else NA_integer_)
      note <- if (is.null(sflist$sf)) "KM unavailable (single/empty group)" else NA_character_

      p_cox <- tryCatch({
        s <- summary(fit$fit)
        as.numeric(s$sctest[3])  # Cox score test p
      }, error = function(e) NA_real_)

      p_mb <- tryCatch({
        postL <- sflist$postL
        id_tab <- table(distinct(postL, id, Z_L)$Z_L)
        two_groups <- length(id_tab) == 2 && all(id_tab > 0)
        has_event  <- sum(postL$event) > 0
        if (!two_groups || !has_event) return(NA_real_)
        cp_for_mb <- postL %>%
          dplyr::transmute(
            tstart, tstop, event,
            response = as.integer(Z_L != ref_label)
          )
        mantel_byar_cp(cp_for_mb, ref_val = 0L)$p
      }, error = function(e) NA_real_)

      tibble(L = fit$L, n_total = fit$n_total, n_events = fit$n_events,
             n_Z0 = z0, n_Z1 = z1,
             HR = fit$hr, CI_low = fit$ci_low, CI_high = fit$ci_high,
             p_cox = p_cox, p_mb  = p_mb,
             note = note)
    }
  }) %>% arrange(L)

  km_plots <- purrr::map2(L_vec, sfs, function(Li, sflist) {
    ttl <- paste0(title_prefix, format(round(Li, 3), trim = TRUE))

    if (is.null(sflist$sf)) {
      return(
        ggplot() +
          annotate("text", x = 0.5, y = 0.5,
                   label = paste0(ttl, "\nKM not available (group size issue)"),
                   size = 4) +
          theme_void() +
          labs(title = if (show_titles) ttl else NULL)
      )
    }

    g <- ggsurvplot(
      sflist$sf,
      data = sflist$postL,
      conf.int = conf_int,
      risk.table = show_risktable,
      risk.table.col = "strata",
      risk.table.y.text = TRUE,
      ggtheme = theme_fn(base_size = base_size),
      palette = palette,
      linetype = "strata",
      censor = TRUE,
      legend.title = NULL,
      legend.labs = c("response=0", "response=1")
    )

    p <- g$plot +
      labs(title = if (show_titles) ttl else NULL,
           x = "Time from L", y = "Survival") +
      theme(legend.position = legend_position,
            plot.title = element_text(face = "bold"),
            legend.title = element_blank()) +
      guides(color   = guide_legend(title = NULL),
             linetype= guide_legend(title = NULL),
             fill    = guide_legend(title = NULL))

    if (!is.null(tmax) && is.finite(tmax)) {
      p <- p + coord_cartesian(xlim = c(0, tmax))
    }

    if (annotate_hrci || annotate_mb_p || annotate_cox_p) {
      row <- tab %>% dplyr::filter(L == Li)
      if (nrow(row) == 1) {
        parts <- character(0)
        HRok <- isTRUE(is.finite(row$HR)) && isTRUE(is.finite(row$CI_low)) && isTRUE(is.finite(row$CI_high))
        if (annotate_hrci && HRok) {
          parts <- c(parts, sprintf("HR=%.2f (%.2f–%.2f)", row$HR, row$CI_low, row$CI_high))
        }
        if (annotate_mb_p && isTRUE(is.finite(row$p_mb))) {
          parts <- c(parts, sprintf("Mantel-Byar test p=%s",
                                    ifelse(row$p_mb < 1e-4, "<1e-4",
                                           formatC(row$p_mb, format = "f", digits = 3))))
        }
        if (annotate_cox_p && isTRUE(is.finite(row$p_cox))) {
          parts <- c(parts, sprintf("Cox p=%s",
                                    ifelse(row$p_cox < 1e-4, "<1e-4",
                                           formatC(row$p_cox, format = "f", digits = 3))))
        }
        if (length(parts) > 0) {
          lab <- paste(parts, collapse = ",\n")
          p <- p + ggplot2::annotate(
            "label", x = Inf, y = Inf,
            label = lab, vjust = 1.1, hjust = 1.1,
            size = 3, label.size = 0.2, lineheight = 0.9
          )
        }
      }
    }

    if (isTRUE(show_risktable) && !is.null(g$table)) {
      tbl <- g$table +
        theme_minimal(base_size = risk_table_base_size) +
        theme(
          legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.margin = margin(t = 0, r = 5, b = 0, l = 5)
        )
      p <- p / tbl + plot_layout(heights = c(1, risk_table_height))
    }
    p
  })

  grid_plot <- wrap_plots(km_plots, ncol = ncol)

  list(
    plot = grid_plot,
    table = tab,
    per_L = list(L = L_vec, fits = fits, survfits = sfs)
  )
}
