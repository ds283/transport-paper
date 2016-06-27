SELECT
  grid.alpha_id                          AS alpha_id,
  grid.beta_id                           AS beta_id,
  grid.alpha                             AS alpha,
  grid.beta                              AS beta,
  grid.tserial                           AS tserial,
  grid.time                              AS time,
  COALESCE(shape_list.shape, -99999999)  AS shape,
  COALESCE(shape_list.redbsp, -99999999) AS redbsp
FROM
  (SELECT
     time_list.tserial AS tserial,
     time_list.time    AS time,
     alpha_list.ROWID  AS alpha_id,
     beta_list.ROWID   AS beta_id,
     alpha_list.alpha  AS alpha,
     beta_list.beta    AS beta
   FROM
     (SELECT DISTINCT (alpha)
      FROM threepf_samples
      ORDER BY alpha) alpha_list,
     (SELECT DISTINCT (beta)
      FROM threepf_samples
      ORDER BY beta) beta_list,
     (SELECT DISTINCT
        (serial) AS tserial,
        time     AS time
      FROM time_samples
      ORDER BY tserial) time_list) grid
  LEFT JOIN
  (SELECT
     threepf_samples.alpha AS alpha,
     threepf_samples.beta  AS beta,
     zeta_threepf.tserial  AS tserial,
     zeta_threepf.threepf  AS shape,
     zeta_threepf.redbsp   AS redbsp
   FROM zeta_threepf
     INNER JOIN threepf_samples ON zeta_threepf.kserial = threepf_samples.serial) shape_list
    ON ABS(shape_list.alpha - grid.alpha) < 1E-4
       AND ABS(shape_list.beta - grid.beta) < 1E-4
       AND shape_list.tserial = grid.tserial
ORDER BY tserial, alpha, beta;
