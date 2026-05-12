# Install required packages if needed
# install.packages(c("hexSticker", "ggplot2", "showtext"))

library(hexSticker)
library(ggplot2)
library(showtext)
library(grid)

# Add a clean font
font_add_google("Montserrat", "montserrat")
font_add_google("Open Sans","open sans")

font_add_google("Roboto","roboto")

navy   <- "#1B3A6B"
blue1  <- "#2166AC"
blue2  <- "#4393C3"
ltblue <- "#AED6F1"
accent <- "#F39C12"   # gold accent for high-precision pan

showtext_auto()

# ============================================================
# BUILD THE BALANCE SCALE PLOT
# ============================================================

create_balance_plot <- function() {

  # --- Fulcrum / pivot point ---
  fulcrum_x <- 0
  fulcrum_y <- 0.10

  # --- Beam: tilted slightly to show "weighting" concept ---
  beam_tilt <- 0.12          # right side slightly lower (higher precision = heavier)
  beam_half  <- 0.75

  beam_left_x  <- fulcrum_x - beam_half
  beam_left_y  <- fulcrum_y + beam_tilt
  beam_right_x <- fulcrum_x + beam_half
  beam_right_y <- fulcrum_y - beam_tilt

  # --- String / suspension from top ---
  string_top_y <- 0.62

  # --- Pan positions (hanging below beam ends) ---
  pan_drop      <- 0.42          # how far pans hang below beam end
  pan_width     <- 0.22
  pan_height    <- 0.04
  pan_curve_pts <- 60

  left_pan_cx  <- beam_left_x
  left_pan_cy  <- beam_left_y  - pan_drop

  right_pan_cx <- beam_right_x
  right_pan_cy <- beam_right_y - pan_drop

  # --- Curve representing a density / precision distribution on each pan ---
  # Left pan: WIDE (low precision, low weight)
  # Right pan: NARROW (high precision, high weight)

  make_curve <- function(cx, cy, sigma, scale_h, n = 200) {
    x_seq <- seq(-0.28, 0.28, length.out = n)
    y_val <- scale_h * exp(-0.5 * (x_seq / sigma)^2)
    data.frame(x = cx + x_seq, y = cy + pan_height / 2 + y_val)
  }

  curve_left  <- make_curve(left_pan_cx,  left_pan_cy,  sigma = 0.18, scale_h = 0.10)
  curve_right <- make_curve(right_pan_cx, right_pan_cy, sigma = 0.07, scale_h = 0.22)

  # --- Shaded area under each curve ---
  shade_left <- rbind(
    data.frame(x = curve_left$x,  y = curve_left$y),
    data.frame(x = rev(curve_left$x),  y = rep(left_pan_cy  + pan_height / 2, nrow(curve_left)))
  )

  shade_right <- rbind(
    data.frame(x = curve_right$x, y = curve_right$y),
    data.frame(x = rev(curve_right$x), y = rep(right_pan_cy + pan_height / 2, nrow(curve_right)))
  )

  # --- Hanging strings from beam ends to pans ---
  string_left  <- data.frame(
    x = c(beam_left_x,  left_pan_cx),
    y = c(beam_left_y,  left_pan_cy + pan_height / 2)
  )
  string_right <- data.frame(
    x = c(beam_right_x, right_pan_cx),
    y = c(beam_right_y, right_pan_cy + pan_height / 2)
  )

  # --- Pan ellipses (draw as thin rectangles with rounded ends via points) ---
  make_pan <- function(cx, cy) {
    theta <- seq(0, pi, length.out = 80)
    data.frame(
      x = cx + (pan_width / 2) * cos(theta),
      y = cy + (pan_height / 2) * sin(theta) + pan_height / 2
    )
  }

  pan_left_arc  <- make_pan(left_pan_cx,  left_pan_cy)
  pan_right_arc <- make_pan(right_pan_cx, right_pan_cy)

  # Pan base lines
  pan_base_left  <- data.frame(
    x = c(left_pan_cx  - pan_width / 2, left_pan_cx  + pan_width / 2),
    y = c(left_pan_cy,  left_pan_cy)
  )
  pan_base_right <- data.frame(
    x = c(right_pan_cx - pan_width / 2, right_pan_cx + pan_width / 2),
    y = c(right_pan_cy, right_pan_cy)
  )

  # --- Fulcrum triangle ---
  tri_h  <- 0.10
  tri_w  <- 0.10
  fulcrum_tri <- data.frame(
    x = c(fulcrum_x - tri_w / 2, fulcrum_x + tri_w / 2, fulcrum_x),
    y = c(fulcrum_y, fulcrum_y, fulcrum_y + tri_h)
  )

  # Stand pole
  pole <- data.frame(
    x = c(fulcrum_x, fulcrum_x),
    y = c(-0.45, fulcrum_y)
  )

  # Stand base
  base_df <- data.frame(
    x = c(-0.20, 0.20),
    y = c(-0.45, -0.45)
  )

  # Suspension string from top of sticker to beam midpoint
  top_string <- data.frame(
    x = c(fulcrum_x, fulcrum_x),
    y = c(string_top_y, fulcrum_y + tri_h)
  )

  # Weight labels
  label_left  <- data.frame(x = left_pan_cx,  y = left_pan_cy  - 0.13,
                            label = "w[i]==frac(1,SE[i]^2)")
  label_right <- data.frame(x = right_pan_cx, y = right_pan_cy - 0.13,
                            label = "w[j]==frac(1,SE[j]^2)")

  # ============================================================
  # ASSEMBLE GGPLOT
  # ============================================================

  navy   <- "#1B3A6B"
  blue1  <- "#2166AC"
  blue2  <- "#4393C3"
  ltblue <- "#AED6F1"
  accent <- "#F39C12"   # gold accent for high-precision pan

  p <- ggplot() +

    # Stand base
    # geom_line(data = base_df, aes(x, y),
    #           color = navy, linewidth = 1.8, lineend = "round") +

    # Pole
    # geom_line(data = pole, aes(x, y),
    #           color = navy, linewidth = 1.4, lineend = "round") +

    # Top suspension string
    geom_line(data = top_string, aes(x, y),
              color = blue1, linewidth = 0.7, linetype = "solid") +

    # Fulcrum triangle
    geom_polygon(data = fulcrum_tri, aes(x, y),
                 fill = blue1, color = blue1, linewidth = 0.5) +

    # Beam
    geom_segment(aes(x = beam_left_x,  y = beam_left_y,
                     xend = beam_right_x, yend = beam_right_y),
                 color = blue1, linewidth = 2.0, lineend = "round") +

    # Hanging strings
    geom_line(data = string_left,  aes(x, y),
              color = blue1, linewidth = 0.8) +
    geom_line(data = string_right, aes(x, y),
              color = blue1, linewidth = 0.8) +

    # Pan arcs (left - low precision)
    geom_path(data = pan_left_arc, aes(x, y),
              color = blue2, linewidth = 1.2) +
    geom_line(data = pan_base_left, aes(x, y),
              color = blue2, linewidth = 1.2) +

    # Pan arcs (right - high precision)
    geom_path(data = pan_right_arc, aes(x, y),
              color = accent, linewidth = 1.4) +
    geom_line(data = pan_base_right, aes(x, y),
              color = accent, linewidth = 1.4) +

    # Shaded density curves
    geom_polygon(data = shade_left, aes(x, y),
                 fill = blue2, alpha = 0.35, color = NA) +
    geom_polygon(data = shade_right, aes(x, y),
                 fill = accent, alpha = 0.40, color = NA) +

    # Density curve lines
    geom_line(data = curve_left, aes(x, y),
              color = blue2, linewidth = 0.9) +
    geom_line(data = curve_right, aes(x, y),
              color = accent, linewidth = 1.1) +

    # Pivot dot
    annotate("point", x = fulcrum_x, y = fulcrum_y + tri_h - 0.1,
             size = 2.5, color = blue1) +

    # Low / High precision labels
    # annotate("text", x = left_pan_cx,  y = left_pan_cy  - 0.06,
    #          label = "Low\nPrecision", size = 2.4, color = blue2,
    #          fontface = "bold", family = "montserrat", lineheight = 0.85) +
    #
    # annotate("text", x = right_pan_cx, y = right_pan_cy - 0.06,
    #          label = "High\nPrecision", size = 2.4, color = accent,
    #          fontface = "bold", family = "montserrat", lineheight = 0.85) +

    coord_fixed(xlim = c(-1.05, 1.05), ylim = c(-0.62, 0.70)) +
    theme_void() +
    theme(
      plot.background  = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA)
    )

  return(p)
}

# ============================================================
# GENERATE THE STICKER
# ============================================================

balance_plot <- create_balance_plot()

sticker(
  subplot      = balance_plot,

  # Package name
  package      = "curveRweights",
  p_size       = 30,
  p_color      = blue1,
  p_family     = "roboto",
  p_fontface   = "plain",
  p_y          = 1.40,           # position text near bottom

  # Subplot placement
  s_x          = 1.00,
  s_y          = 0.85,
  s_width      = 1.55,
  s_height     = 1.35,

  # Hex colours
  h_fill       = "#FFFFFF",      # white background
  h_color      = blue1,      # navy blue border
  h_size       = 2.0,

  # Output
  filename     = "man/figures/logo.png",
  dpi          = 600,
  white_around_sticker = FALSE
)

sticker(
  subplot      = balance_plot,

  # Package name
  package      = "curveRweights",
  p_size       = 8,
  p_color      = blue1,
  p_family     = "roboto",
  p_fontface   = "plain",
  p_y          = 1.40,           # position text near bottom

  # Subplot placement
  s_x          = 1.00,
  s_y          = 0.85,
  s_width      = 1.55,
  s_height     = 1.35,

  # Hex colours
  h_fill       = "#FFFFFF",      # white background
  h_color      = blue1,      # navy blue border
  h_size       = 2.0,

  # Output
  filename     = "man/figures/logo_small.png",
  dpi          = 150,
  white_around_sticker = FALSE
)

message("✅  Hex sticker saved as  cman/figures/logop.ng")

