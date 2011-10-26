DROP TABLE IF EXISTS `map`;
CREATE TABLE `map` (
  `id` MEDIUMINT UNSIGNED NOT NULL PRIMARY KEY AUTO_INCREMENT,
  `name` VARCHAR(255),
  `xsize` TINYINT UNSIGNED,
  `ysize` TINYINT UNSIGNED,
  `bootstrap` SMALLINT UNSIGNED,
  `commandline` TEXT,
  UNIQUE KEY(`name`, `bootstrap`)
);

DROP TABLE IF EXISTS `vecdim`;
CREATE TABLE `vecdim` (
  `id` INT UNSIGNED NOT NULL PRIMARY KEY AUTO_INCREMENT,
  `map` MEDIUMINT UNSIGNED NOT NULL,
  `idx` SMALLINT UNSIGNED,
  `count_pos` SMALLINT UNSIGNED,
  `name` TEXT,
   KEY(`map`, `idx`)
);

DROP TABLE IF EXISTS `mapnode`;
CREATE TABLE `mapnode` (
  `id` INT UNSIGNED NOT NULL PRIMARY KEY AUTO_INCREMENT,
  `map` MEDIUMINT UNSIGNED NOT NULL,
  `x` TINYINT UNSIGNED NOT NULL,
  `y` TINYINT UNSIGNED NOT NULL,
   KEY(`map`, `x`, `y`)
);

DROP TABLE IF EXISTS `weight`;
CREATE TABLE `weight` (
  `id` INT UNSIGNED NOT NULL PRIMARY KEY AUTO_INCREMENT,
  `mapnode` INT UNSIGNED NOT NULL,
  `vecdim` INT UNSIGNED NOT NULL,
  `weight` FLOAT,
  KEY(`vecdim`),
  KEY(`mapnode`, `vecdim`)
);
  

DROP TABLE IF EXISTS `mapping`;
CREATE TABLE `mapping` (
  `id` INT UNSIGNED NOT NULL PRIMARY KEY AUTO_INCREMENT,
  `map` MEDIUMINT UNSIGNED NOT NULL,
  `mapnode` INT UNSIGNED NOT NULL,
  `gene_id` VARCHAR(50) NOT NULL,
  `error` FLOAT NOT NULL,
  KEY(`map`, `mapnode`),
  KEY(`mapnode`),
);

DROP TABLE IF EXISTS `bit`;
CREATE TABLE `bit` (
  `id` INT UNSIGNED NOT NULL PRIMARY KEY AUTO_INCREMENT,
  `mapping` INT UNSIGNED NOT NULL,
  `vecdim` INT UNSIGNED NOT NULL,
  `state` enum('nan', '0', '1'),
  KEY(`mapping`, `vecdim`)
);

DROP TABLE IF EXISTS `expression`;
CREATE TABLE `expression` (
  `id` INT UNSIGNED NOT NULL PRIMARY KEY AUTO_INCREMENT,
  `mapping` INT UNSIGNED NOT NULL,
  `vecdim` INT UNSIGNED NOT NULL,
  `value` FLOAT,
  KEY(`mapping`, `vecdim`)
);

DROP TABLE IF EXISTS `go_ora_result`;
CREATE TABLE `go_ora_result` (
  `id` INT UNSIGNED NOT NULL PRIMARY KEY AUTO_INCREMENT,
  `mapnode` INT UNSIGNED NOT NULL,
  `term` varchar(255) NOT NULL,
  `score` FLOAT NOT NULL,
  `numgenes` MEDIUMINT UNSIGNED NOT NULL,
  `corrected_pvalue` FLOAT NOT NULL,
  `pvalue` FLOAT NOT NULL,
  `serial_number` SMALLINT UNSIGNED NOT NULL,
  KEY(`mapnode`)
);
  

