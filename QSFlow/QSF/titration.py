import subprocess
import os

class titration:
    def __init__(self, titration_matrix, key, outputdir, solvent, solute,path_inital, intial=False) :
        titration_list = []
        titration_list[:] = titration_matrix
        titration_list.pop(titration_matrix.index(1.0))
        print(intial)


        if not intial:
            for i in range(len(titration_list)):

                subprocess.run(
                    [f"cp -r {os.path.join(outputdir,f'InputGrofiles{key}')} {os.path.join(outputdir, f'InputGrofiles{key}_{titration_list[i]}')}"],
                    shell=True)

                for iteams in solute:
                    chargeMatrix = []
                    with open(
                            f'{os.path.join(outputdir, f"InputGrofiles{key}_{titration_list[i]}", f"{iteams[:3]}_Solute1.itp")}', 'r') as gmx:
                        lines = gmx.readlines()

                        index_GMXESP = 0
                        index_lastGMXESP = 0
                        for line in lines:
                            if line.strip() == "[ atoms ]":
                                index_GMXESP = lines.index(line) + 2
                                break
                        for line in lines:
                            if line.strip() == "[ bonds ]":
                                index_lastGMXESP = lines.index(line) - 2

                        for a in range(index_lastGMXESP - index_GMXESP + 1):
                            charge_int = float(lines[index_GMXESP + a].split()[6]) * float(titration_list[i])
                            charge = f'{charge_int:.4f}'.rjust(11)

                            chargeMatrix.append(charge)
                        for char in range(index_lastGMXESP - index_GMXESP + 1):
                            chargeIndexGmx = lines[index_GMXESP - 1].index("charge") - 6

                            lines[index_GMXESP + char] = lines[index_GMXESP + char][0:chargeIndexGmx] + str(
                                chargeMatrix[char]) + lines[index_GMXESP + char][chargeIndexGmx + 11:]
                        subprocess.run(
                            [f'rm {os.path.join(outputdir, f"InputGrofiles{key}_{titration_list[i]}", f"{iteams[:3]}_Solute1.itp")}'],
                            shell=True)

                    with open(f'{os.path.join(outputdir, f"InputGrofiles{key}_{titration_list[i]}", f"{iteams[:3]}_Solute1.itp")}',
                              'a') as itp:
                        for line in lines:
                            itp.writelines(line)

        else:
            for i in range(len(titration_list)):
                subprocess.run(
                    [
                        f"cp -r {os.path.join(outputdir, f'InputGrofiles{key}')} {os.path.join(outputdir, f'InputGrofiles{key}_{titration_list[i]}')}"],
                    shell=True)

                for iteams in solute:
                    chargeMatrix = []
                    with open(
                            f'{os.path.join(outputdir, f"InputGrofiles{key}_{titration_list[i]}", f"{iteams[:3]}_Solute1.itp")}',
                            'r') as gmx:
                        lines = gmx.readlines()

                        index_GMXESP = 0
                        index_lastGMXESP = 0
                        for line in lines:
                            if line.strip() == "[ atoms ]":
                                index_GMXESP = lines.index(line) + 2
                                break
                        for line in lines:
                            if line.strip() == "[ bonds ]":
                                index_lastGMXESP = lines.index(line) - 2

                        for a in range(index_lastGMXESP - index_GMXESP + 1):
                            charge_int = float(lines[index_GMXESP + a].split()[6]) * float(titration_list[i])
                            charge = f'{charge_int:.4f}'.rjust(11)

                            chargeMatrix.append(charge)
                        for char in range(index_lastGMXESP - index_GMXESP + 1):
                            chargeIndexGmx = lines[index_GMXESP - 1].index("charge") - 6

                            lines[index_GMXESP + char] = lines[index_GMXESP + char][0:chargeIndexGmx] + str(
                                chargeMatrix[char]) + lines[index_GMXESP + char][chargeIndexGmx + 11:]
                        subprocess.run(
                            [
                                f'rm {os.path.join(outputdir, f"InputGrofiles{key}_{titration_list[i]}", f"{iteams[:3]}_Solute1.itp")}'],
                            shell=True)

                    with open(
                            f'{os.path.join(outputdir, f"InputGrofiles{key}_{titration_list[i]}", f"{iteams[:3]}_Solute1.itp")}',
                            'a') as itp:
                        for line in lines:
                            itp.writelines(line)