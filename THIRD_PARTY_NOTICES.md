# Third-party notices

## mpmath: de Hoog quotient-difference recurrence

`matlab/private/response_dehoog.m` adapts the quotient-difference table,
continued-fraction recurrence and improved remainder from the `deHoog`
class in mpmath's `mpmath/calculus/inverselaplace.py`, consulted 2026-09-26:
https://github.com/mpmath/mpmath/blob/master/mpmath/calculus/inverselaplace.py

That source credits Kristopher L. Kuhlman (February 2017). The MATLAB
domain policy, fixed-double-precision safeguards, adaptive comparisons,
diagnostics and API are separate additions. No Python runtime is required.
The recurrence implements de Hoog, Knight and Stokes (1982),
https://doi.org/10.1137/0903022.

The adapted code is distributed under the following BSD-3-Clause terms:

Copyright (c) 2005-2026 Fredrik Johansson and mpmath contributors

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
  a. Redistributions of source code must retain the above copyright notice,
     this list of conditions and the following disclaimer.
  b. Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the following disclaimer in the
     documentation and/or other materials provided with the distribution.
  c. Neither the name of the copyright holder nor the names of its
     contributors may be used to endorse or promote products derived
     from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
