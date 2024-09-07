import subprocess
import os


class toopol:
    """
    A class to create topol.top and nmol.top files for the MD calculation.

    """

    def __init__(
        self,
        solvent1=None,
        solvent1N=None,
        solvent2="",
        solvent2N=None,
        solute=None,
        con=None,
        currentdir=None,
        x=None,
        y=None,
        z=None,
        key=None,
        inital=False,
    ):
        """

        :param solvent1: Name of solvent 1.
        :param solvent1N: Name of solvent 1.
        :param solvent2: NOT IN USE, will be used in the future.
        Name of the second solvent for multi solvent simulations.
        :param solvent2N: NOT IN USE, will be used in the future.
        Number of the second solvent for multi solvent simulations.
        :param solute: Solute matrix.
        :param con: Solute concentration matrix.
        :param currentdir: Current directory of input files.
        :param x: X-length of topology.
        :param y:  Y-length of topology.
        :param z:  Z-length of topology.
        :param key: Key to the current system.
        :param inital: If initial system is given or not.


        """
        self.solvent = solvent1
        self.solvent2 = solvent2
        self.solute1 = solute
        self.con = con
        self.dir = currentdir
        self.x = float(x) if x is not None else None
        print(self.x)
        self.y = float(y) if y is not None else None
        self.z = float(z) if z is not None else None
        command2 = f'mkdir {os.path.join(self.dir, f"InputGrofiles{key}")}'
        if os.path.isdir(os.path.join(self.dir, f"InputGrofiles{key}")):
            folder_number = 0
            if not os.path.isdir(os.path.join(self.dir, "Old_InputGrofiles_Output")):
                make_old = f'mkdir {os.path.join(self.dir, "Old_InputGrofiles_Output")}'
                subprocess.run(make_old, shell=True, check=True)
            while True:
                if os.path.isdir(
                    os.path.join(
                        self.dir,
                        f"Old_InputGrofiles_Output",
                        f"InputGrofiles{folder_number + 1}{key}",
                    )
                ):
                    folder_number += 1
                else:
                    move = f'mv {os.path.join(self.dir, f"InputGrofiles{key}")} {os.path.join(self.dir, f"Old_InputGrofiles_Output", f"InputGrofiles{folder_number + 1}{key}")}'
                    subprocess.run(move, shell=True, check=True)
                    break
        if os.path.isdir(os.path.join(self.dir, f"Output{key}")):
            folder_number = 0
            while True:
                if os.path.isdir(
                    os.path.join(
                        self.dir,
                        f"Old_InputGrofiles_Output",
                        f"Output{folder_number + 1}{key}",
                    )
                ):
                    folder_number += 1
                else:
                    move = f'mv {os.path.join(self.dir, f"Output{key}")} {os.path.join(self.dir, f"Old_InputGrofiles_Output", f"Output{folder_number + 1}{key}")}'
                    subprocess.run(move, shell=True, check=True)
                    break
        subprocess.run(command2, shell=True, check=True)
        if not inital:
            cmd = f"touch {os.path.join(self.dir, f'InputGrofiles{key}', 'topol.top')} && touch {os.path.join(self.dir, f'InputGrofiles{key}', 'nmol.itp')}"
            subprocess.run(cmd, shell=True, check=True)
            with open(
                os.path.join(self.dir, f"InputGrofiles{key}", "topol.top"), "a"
            ) as top, open(
                os.path.join(self.dir, f"InputGrofiles{key}", "nmol.itp"), "a"
            ) as nmol:

                lines_atomtypes = [
                    '#include "oplsaa.ff/forcefield.itp"',
                    f'#include "{self.solvent}_Solvent_atomtype.itp"',
                ]
                lines_itp = [f'#include "{self.solvent}_Solvent.itp"']

                last_lines = ["", "[system]", f"{solvent1}", "", '#include "nmol.itp"']

                if len(self.solvent2) != 0:
                    lines_atomtypes.append(
                        f'#include "{self.solvent2}_Solvent2_atomtype.itp"'
                    )
                    lines_itp.append(f'#include "{self.solvent2}_Solvent2.itp"')
                for i in range(len(self.solute1)):
                    lines_atomtypes.append(
                        f'#include "{self.solute1[i].strip()[:3]}_Solute{1}_atomtype.itp"'
                    )
                    lines_itp.append(
                        f'#include "{self.solute1[i].strip()[:3]}_Solute{1}.itp"'
                    )

                for l in lines_atomtypes:
                    top.write(l + "\n")
                for sha in lines_itp:
                    top.write(sha + "\n")
                for lami in last_lines:
                    top.write(lami + "\n")

                linesNmol = [
                    "[ molecules ]",
                    ";molecules #molecules",
                    "",
                    f"{self.solvent}   {solvent1N}",
                ]
                if len(self.solvent2) != 0:
                    linesNmol.append(f"{self.solvent2}   {solvent2N}")

                for i in range(len(self.solute1)):
                    self.mol = int(
                        float(self.con[i].strip())
                        * 6.4e-20
                        * 6.02214e23
                        * ((self.x * self.y * self.z) / (40 * 40 * 40))
                    )
                    linesNmol.append(f"{self.solute1[i].strip()[:3]}   {self.mol}")

                for li in linesNmol:
                    nmol.write(li + "\n")
