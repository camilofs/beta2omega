# beta2omega
A script to convert beta (bcc) to omega (hcp), common crystal structures among Ti alloys. <p>

The code generates four distinct omega-phase structures based on a beta (bcc) input. The input must be in the VASP (POSCAR) format. The functions are name as beta_to_omegaN, where N = (1..4). Usage:

```python
import beta2omega
from beta2omega import *

_header, _atoms = load_atoms()
_omega4 = beta_to_omega4(_atoms) # fourth variant
print_to_file(_header, _omega4)
```

### Details on the beta -> transformation
The script works with a 3x3x3 bcc superell and direct coordinates (0-1); in this case, the smallest coord. would be 1/6
i.e. the position of the first centered atom in the bcc structure is at (1/6, 1/6, 1/6). For simplification purposes, the code treats coordinates as a f(a) = 1/a (thus the first atom is at 6,6,6).

To visualize the operations coded in this script, open the POSCAR file in VESTA - **where [p] = projection vector and [u] upward vector** - and follow the instructions below:

### Omega variant 1
In this variant the sums of coordinates determine the shift **@ref1**. <p>
(VESTA alignment [p]1-10 [u]111, at @ 0,0,0 moves upwards):

    Invar       02,05,08,11,14
    Up  -->   0,03,06,09,12,15,18
    Down <-     04,07,10,13,16

### Omega variant 2
Similar to variant 1. <p>
(VESTA alignment [p]1-10 [u]111, at @ 0,0,0 moves down instead):

    Invar    04,07,10,13,16
    Up -->   02,05,08,11,14
    Down - 0,03,06,09,12,15,18


### Omega variant 3
In this variant, the [z] axis ascertain the shift, with a few specific conditions **@ref2**. <p>
VESTA alignment [p]110 [u]1-11, atom @ 0,1,0 moves upwards

**Condition 1: a = b = c**

    Invar       02,05
    Up  -->  00,03,06
    Down <-     01,04

Due to PBC, 0 == 6 and; if a = b

    Invar       ab2,ab5
    Up  --> ab0,ab3,ab6
    Down <-     ab1,ab4

**Condition 2: else:**
Due to PBC (0-6), 2i = 4
The remaining condition is:

For z = 2, 5

    Down <- if a-b == 2
    Up  --> if a-b == 4

For z = 0, 3, 6

    Invar   if a-b == 2
    Down <- if a-b == 4

For z = 1, 4

    Up -->  if a-b == 2
    Invar   if a-b == 4

**@ref3**
A final adjustment is needed when aplying the shift to the structure.
Since the upward vector is <1-11>, we must apply an (a-bc) shift to preserve our convention of a positive z axis.

### Omega variant 4 
Similar to variant #3. <p>
VESTA alignment [p]110 [u]1-11, atom @ 0,1,0 moves down instead, thus:

    Invar       01,04
    Up  -->     02,05
    Down <-  00,03,06

and so on (...)

