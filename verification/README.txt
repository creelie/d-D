Ancillary files for "Class Numbers up to Ten Thousand and Zero-Free Half-Planes" (Deep Bhattacharjee).
lbound.py, uniform.py, opt.py, opt2.py: certified bounds (python-flint / Arb, 128 bits).
  python3 uniform.py   -> constants of Theorem 5.2 (Appendix B)
classcount.c: exact class numbers by counting reduced forms.
  gcc -O3 -fopenmp -o classcount classcount.c -lm
  OMP_NUM_THREADS=4 ./classcount 339110855 1000 classnumbers_1000.txt   (about 15 minutes)
classnumbers_1000.txt: header line, then "h F(h) max|d|" for h = 1..1000.
hcheck.py: independent check of single class numbers.
mkdata.py, fig*.py: tables and figures of the paper.
sanity.py: numerical check that the weighted prime sum of Corollary 2.3 matches log L(1,chi_d).
tatuzawa.py: the unconditional bound h(d) >= 2901 of Theorem 1.4 (one possible exception).
