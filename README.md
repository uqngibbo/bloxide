# bloxide -- A compressible boundary layer analysis code

bloxide is a tool for analysing self-similar compressible boundary layers, written in Rust.

## Build Instructions
Install dependancies:

    curl https://sh.rustup.rs -sSf | sh

Clone the repository:

     git clone https://github.com/uqngibbo/bloxide.git

Build the code:

    cd bloxide
    cargo install --path .

## Example Use

    cd examples
    bloxide cold.yaml

This should output some details about the boundary layer as well as a file called cold.dat which has the entire profile saved in it.

    bloxide: A compressible boundary layer analysis code.
    Config {
        R: 287.1,
        gamma: 1.4,
        Pr: 0.71,
        p_e: 2303.0,
        u_e: 604.5,
        T_e: 108.1,
        T_wall: 269.5,
        x: 0.5,
    }
    Solved fdd 0.5013822034038323 gd -0.0363666169277437 in 5 iters
    Skin Friction:  5.05226 N/m2
    Heat Transfer : -0.00927 W/cm2
    99.9% BL size : 2.77701 mm
    Rex: 2.99296 million   Ret: 16622.97552
    Solved fdd 0.4988390188851437 g 2.411558697875749 in 5 iters
    Adiabatic Wall Temp: 260.68950 K
    Done.

## Author:
Nick Gibbons (n.gibbons@uq.edu.au) and Peter Jacobs

## License:
bloxide use is governed by the GNU General Public License 3. This is a copyleft license, which means that any code based on bloxide must also be made freely available under a similar license.
See the file gpl.txt for further details.
