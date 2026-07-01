# Color utilities

Bulletproof color handling for figure generation. Validates hex codes,
falls back gracefully on invalid inputs, and builds n-color vectors that
never crash downstream rendering (col2rgb / scale_color_manual).
