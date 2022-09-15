import pandas as pd
import re
import sympy as sm

# Import data about elements
interests = ['Symbol', 'Element', 'AtomicNumber', 'AtomicMass']
elements = pd.read_csv('periodic_table.csv')[interests]


def gcd(a, b):
    if (b == 0):
        return a
    else:
        return gcd(b, a % b)


def lcm(a, b): return (a*b)//gcd(a, b)


def listlcm(L):
    res = L[0]
    for n in L[1::]:
        res = lcm(n, res)
    return res


def literalize(eqn):
    result = eqn
    # Convert all parenthesis multiples into element forms. For example, (CH2)6 becomes C6H12.
    # Maintain groupings with the ~ tag (such as polyatomic ions that must remain intact throughout the reaction process)
    step = re.findall("\(((?![~])[A-Za-z0-9]+)\)([0-9]+)", eqn)
    repls = []
    for match in step:
        rep = ""
        split1 = re.findall("[A-Z][a-z]?[0-9]*", match[0])
        for elem in split1:
            symb = re.findall("[A-Z][a-z]?", elem)[0]
            count = re.findall("\d+", elem)
            if len(count) < 1:
                count.append("1")
            num = int(count[0]) * int(match[1])
            rep += symb + str(num)
        result = re.sub(f"\(~*({match[0]})\)({match[1]})", f"{rep}", result)
    return result


def balance(eqn, lit=True, star=False, unk=False, el=False, mat=False, ones=False):
    if lit:
        eqn = literalize(eqn)
    lhs, rhs = re.split("=", eqn)
    lhl = re.split("[+]", lhs)
    rhl = re.split("[+]", rhs)
    elems = []
    for i in range(len(lhl)):
        lhl[i] = re.findall("[A-Z][a-z]?[0-9]*|\(~\w+\)[0-9]*", lhl[i])
        for j in range(len(lhl[i])):
            symb = re.findall("[A-Z][a-z]?|\(~\w+\)", lhl[i][j])[0]
            if not symb in elems:
                elems.append(symb)
            count = re.findall("[0-9]+", lhl[i][j])
            if len(count) < 1:
                count.append("1")
            lhl[i][j] = [symb, int(count[0])]
    for i in range(len(rhl)):
        rhl[i] = re.findall("[A-Z][a-z]?[0-9]*|\(~\w+\)[0-9]*", rhl[i])
        for j in range(len(rhl[i])):
            symb = re.findall("[A-Z][a-z]?|\(~\w+\)", rhl[i][j])[0]
            if not symb in elems:
                elems.append(symb)
            count = re.findall("[0-9]+", rhl[i][j])
            if len(count) < 1:
                count.append("1")
            rhl[i][j] = [symb, int(count[0])]
    if el:
        print(elems)
    matrix = []
    for i in range(len(elems)):
        current = elems[i]
        inmat = []
        for j in range(len(lhl)):
            entered = False
            for k in range(len(lhl[j])):
                if lhl[j][k][0] == current:
                    if entered:
                        inmat[-1] += lhl[j][k][1]
                    else:
                        inmat.append(lhl[j][k][1])
                        entered = True
            if not entered:
                inmat.append(0)
        for j in range(len(rhl)):
            entered = False
            for k in range(len(rhl[j])):
                if rhl[j][k][0] == current:
                    if entered:
                        inmat[-1] -= rhl[j][k][1]
                    else:
                        inmat.append(-rhl[j][k][1])
                        entered = True
            if not entered:
                inmat.append(0)
        matrix.append(inmat)
    if mat: print(matrix)
    reduced = (sm.Matrix(matrix)).rref()
    if mat: print(reduced)
    endcol = reduced[0].col(-1)
    denoms = []
    for i in range(len(endcol)):
        denoms.append(endcol[i].as_numer_denom()[1])
    factor = listlcm(denoms)
    unknowns = []
    for i in range(len(lhl)+len(rhl)-1):
        unknowns.append(-factor*endcol[i])
    unknowns.append(factor)
    if unk:
        print(unknowns)
    outstr = ""
    for i in range(len(lhl)):
        if ((unknowns[i] != 1) or (ones == True)):
            outstr += str(unknowns[i])
            if star:
                outstr += "*"
        for j in range(len(lhl[i])):
            outstr += lhl[i][j][0]
            if lhl[i][j][1] > 1:
                outstr += str(lhl[i][j][1])
        if i == len(lhl) - 1:
            outstr += "="
        else:
            outstr += "+"
    for i in range(len(rhl)):
        if ((unknowns[len(lhl)+i] != 1) or (ones == True)):
            outstr += str(unknowns[len(lhl)+i])
            if star:
                outstr += "*"
        for j in range(len(rhl[i])):
            outstr += rhl[i][j][0]
            if rhl[i][j][1] > 1:
                outstr += str(rhl[i][j][1])
        if not i == len(rhl) - 1:
            outstr += "+"
    return outstr


def details(eqn):  # Input a balanced equation, and this will describe some parts of the reaction
    print("Details on reaction: " + eqn)
    lhs, rhs = re.split("=", eqn)
    lhl = re.split("[+]", lhs)
    rhl = re.split("[+]", rhs)
    out1 = "Reactants: "
    out2 = "Products: "
    out3 = "Reactant mass composition: "
    out4 = "Product mass composition: "
    masses = []
    try:
        for i in range(len(lhl)):  # Look at each compound on the LHS
            split = re.findall("([0-9]|^)([A-Za-z0-9]+)", lhl[i])[0]
            coeff = split[0]
            if coeff == "":
                coeff = "1"
            coeff = int(coeff)
            mass = 0
            elsplit = re.findall("([A-Z][a-z]?)([0-9]*)", split[1])
            for elem in elsplit:
                count = elem[1]
                if count == "":
                    count = "1"
                count = int(count)
                mass += float(elements[elements.Symbol ==
                              elem[0]].AtomicMass) * count * coeff
            masses.append(mass)
            out1 = f"{out1} {str(split[0])}N mol [{mass:.3f}N g] {split[1]}"
            if i != len(lhl) - 1:
                out1 += ", "
        for i in range(len(rhl)):  # Look at each compound on the RHS
            split = re.findall("([0-9]+|^)([A-Za-z0-9]+)", rhl[i])[0]
            coeff = split[0]
            if coeff == "":
                coeff = "1"
            coeff = int(coeff)
            mass = 0
            elsplit = re.findall("([A-Z][a-z]?)([0-9]*)", split[1])
            for elem in elsplit:
                count = elem[1]
                if count == "":
                    count = "1"
                count = int(count)
                mass += float(elements[elements.Symbol ==
                              elem[0]].AtomicMass) * count * coeff
            masses.append(mass)
            out2 = f"{out2} {str(split[0])}N mol [{mass:.3f}N g] {split[1]}"
            if i != len(rhl) - 1:
                out2 += ", "
        totalmass = 0
        for mass in masses:
            totalmass += mass
        totalmass /= 2
        for i in range(len(lhl)):
            out3 = f"{out3}{100*masses[i]/totalmass:.3f}%"
            if i != len(lhl) - 1:
                out3 += ", "
        for i in range(len(rhl)):
            out4 = f"{out4}{100*masses[len(lhl)+i]/totalmass:.3f}%"
            if i != len(rhl) - 1:
                out4 += ", "
    except KeyError:
        print("Could not find a referenced element.")
    except Exception as e:
        print("An error occurred in determining details for the given reaction.")
        print(e)
    else:
        print(out1)
        print(out2)
        print(out3)
        print(out4)
    return

def composition(comp, printmass=True): # Takes a formula for a chemical compound and prints the percent composition, by mass, of its constituent elements
    elems = [[],[]]
    totalmass = 0
    outstr = f"{comp} elemental mass compositions: "
    out2 = ""
    try:
        elsplit = re.findall("([A-Z][a-z]?)([0-9]*)",comp)
        for elem in elsplit:
            name = elem[0]
            count = elem[1]
            if count == "": count = "1"
            x = float(elements[elements.Symbol == name].AtomicMass) * int(count)
            totalmass += x
            if name in elems[0]:
                i = elems[0].index(name)
                elems[1][i] += x
            else:
                elems[0].append(name)
                elems[1].append(x)
        for i in range(len(elems[0])):
            if i != 0: outstr += ", "
            outstr += f"{100*elems[1][i]/totalmass:.3f}% {elems[0][i]}"
        out2 = f"  1 mol has an atomic mass of {totalmass:.3f} g"
    except KeyError:
        print("Could not find a referenced element.")
    except Exception as e:
        print("An error occurred in calculating the elemental mass compositions of the given compound.")
        print(e)
    else:
        print(outstr)
        if printmass: print(out2)
    return
