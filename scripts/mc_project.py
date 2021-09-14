import argparse
import random
import math
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Programme permettant de realiser des repliement proteique en se basant'
                                             'sur une methode de Monte Carlo avec un calcul energetique '
                                             'sur les acides amines hydrophobes ')
parser.add_argument('sequence', nargs=1, metavar='sequence', type=str, help="The sequence to be used.")
parser.add_argument('-t', '--temperature', nargs=1, help="choose the temperature for the Monte-Carlo search"
                                                         "(temperature advised between 0.01 and 0.9)")
parser.add_argument('-s', '--step', nargs=1, help="nombre d'etape du MCsearch (par defaut 500)")

args = parser.parse_args()

AAapolaire = ["A", "I", "L", "M", "F", "P", "V"]
AApolaire = ["D", "E", "R", "N", "C", "Q", "G", "H", "K", "S", "T", "W", "Y"]


class AminoAcid(object):
    envalue = 0
    coor_i = 0
    coor_j = 0

    def __init__(self, name, pol, pred, succ, hid_name):
        self.name = name
        self.pol = pol
        self.pred = pred
        self.succ = succ
        self.hid_name = hid_name

    def get_name(self):
        return self.name

    def get_hid_name(self):
        return self.hid_name

    def get_pol(self):
        return self.pol

    def get_pred(self):
        return self.pred

    def get_succ(self):
        return self.succ

    def set_value(self, value):
        self.envalue = value

    def get_value(self):
        return self.envalue

    def set_pred(self, pred):
        self.pred = pred

    def update_coor(self, i, j):
        self.coor_i = i
        self.coor_j = j

    def get_coor(self):
        return self.coor_i, self.coor_j

    def show(self):
        print(f"pol: {self.pol}, pred: {self.pred}, succ: {self.succ}")

    def show_coor(self):
        print(f"{self.name} : i = {self.coor_i}  j = {self.coor_j}")


def translation(seq, apol):
    """
    Permet de traduire une sequence selon la polarite de ces acides amines en transformant ces acides amines en objet
    AminoAcid et en conservant ainsi sa polarite et sa position dans la chaine proteique
    :param seq: sequence d'interet
    :param apol: liste des acide amine hydrophobe
    :return: un dictionnaire d'objet d'AminoAcid avec pour cle la position de l'acide amine dans la sequence
    """
    trans_seq = {}
    count = 0
    for aa in seq:
        if aa in apol:
            if count != 0:
                trans_seq[count] = AminoAcid(count, "H", count - 1, count + 1, aa)
            else:
                trans_seq[count] = AminoAcid(count, "H", "none", count + 1, aa)
        else:
            if count != 0:
                trans_seq[count] = AminoAcid(count, "P", count - 1, count + 1, aa)
            else:
                trans_seq[count] = AminoAcid(count, "P", "none", count + 1, aa)
        count += 1
    return trans_seq


def mat_maker(size):
    """
    permet de creer une matrice vide (rempli de 0)
    :param size: taille de la matrice
    :return: une matrice
    """
    mat = [[0] * size for i in range(size)]
    return mat


def mat_filler(mat, seq):
    """
    Permet de remplir la matrice avec les AminoAcid de la sequence aligne au centre de la matrice
    :param mat: matrice vide (rempli de 0)
    :param seq: dictionnaire d'AminoAcid dans le bon ordre de la sequence
    :return: une matrice avec la sequence d'interet a son centre
    """
    mid = len(mat) / 2
    mid = int(mid)
    half_seq = len(seq) / 2
    half_seq = int(half_seq)
    for aa in seq:
        mat[mid][mid - half_seq] = seq[aa]
        seq[aa].update_coor(mid, mid - half_seq)
        half_seq -= 1
    return mat


def printer(mat):
    """
    Fonction permettant d'afficher une matrice en ASCII en parcourant la dite matrice
    :param mat: une matrice
    :return: ne renvoie rien mais imprime
    """
    for i in mat:
        for j in i:
            if type(j) == int:
                print(f"  {j} ", end='')
            else:
                to_print = str(j.get_name()) + j.get_pol()
                print(f" {to_print} ", end='')
        print("")


def check_neighbors(mat, aa, bina):
    """
    Fonction permettant de trouver les voisins vertical et horizontal d'un acide amine et de conserver les positions
    des 0 dans une liste et les acides amines dans une autre et de renvoyer l'une des deux listes selon la valeur du
    parametre bina
    :param mat: la matrice a explorer
    :param aa: AminoAcid dont on recherche les voisins dans la matrice
    :param bina: 1 si on veut la liste des AminoAcid en sortie et 0(ou autre valeur que 1) si on veut la liste des 0
    :return: une liste d'AminoAcid ou une liste de position de 0
    """
    i, j = aa.get_coor()
    neigh = []
    neighnull = []
    turn = [-1, 1]
    for posturn in turn:
        if mat[i + posturn][j] != 0:
            neigh.append(mat[i + posturn][j])
        else:
            neighnull.append((i + posturn, j))
        if mat[i][j + posturn] != 0:
            neigh.append(mat[i][j + posturn])
        else:
            neighnull.append((i, j + posturn))
    if bina == 1:
        return neigh
    else:
        return neighnull


def vshd_move_ends(neighbor, mat, aa):
    """
    Permet d'effectuer un movement VSHD pour les AminoAcid aux extremites de la sequence en trouvant les positions
    voisines ou un deplacement est possible, si plusieurs deplacement sont possible un de ces deplacements est choisi
    au hasard
    :param neighbor: liste de voisin
    :param mat: matrice avec des AminoAcid
    :param aa: AminoAcid d'interet
    :return: ne renvoie rien modifie directement la matrice ainsi que les coordonnees de l'objet
    """
    aai, aaj = aa.get_coor()
    empty_slots = check_neighbors(mat, neighbor, 0)
    if len(empty_slots) == 0:
        return 0
    choice = random.choice(empty_slots)
    mat[choice[0]][choice[1]] = aa
    aa.update_coor(choice[0], choice[1])
    mat[aai][aaj] = 0


def direct_neighbors(neighbors, target):
    """
    Permet a partir d'une liste de voisin d'un AminoAcid de determiner les voisins qui lui sont directement liee dans
    la sequence en utilisant les propriete des objets AminoAcid (leurs noms correspondant a leurs positions dans la
    sequence)
    :param neighbors: liste de voisin d'un AminoAcid d'interet
    :param target: AminoAcid d'interet
    :return: une liste d'AminoAcid (2) directement lie a l'AminoAcid d'interet
    """
    linked_neighbors = []
    for neighbor in neighbors:
        neiname = neighbor.get_name()
        if neiname == target + 1 or neiname == target - 1:
            linked_neighbors.append(neighbor)
    return linked_neighbors


def vshd_move_base(neighbors, mat, aa):
    """
    Permet d'effectuer un deplacement VSHD en diagonale d'un AminoAcid lorsque ce deplacement est possible en verifiant
    par rapport a ses AminoAcid lie si cela est possible
    :param neighbors: liste de voisins de l'AminoAcid d'interet
    :param mat: matrice d'AminoAcid
    :param aa: AminoAcid d'interet
    :return: 0 et modifie directement la matrice et les coordonnees de l'objet AminoAcid si le deplacement est possible
    """
    aai, aaj = aa.get_coor()
    check = 0
    aaname = aa.get_name()
    linked_neighbors = direct_neighbors(neighbors, aaname)
    empty1 = check_neighbors(mat, linked_neighbors[0], 0)
    empty2 = check_neighbors(mat, linked_neighbors[1], 0)
    for slot in empty1:
        if slot in empty2:
            mat[slot[0]][slot[1]] = aa
            aa.update_coor(slot[0], slot[1])
            mat[aai][aaj] = 0
            check = 1
    if check == 0:
        return False
    else:
        return True


def vshd_crankshaft_test(mat, aa, prots):
    """
    Cette fonction permet de verifier si un deplacement de crankshaft est possible pour un AminoAcid en verifiant que
    toutes les conditions sont bien reuni comme le fait qu'on ait bien une formation d'AminoAcid en 'carre' et deux
    espaces vides (0) en face des AminoAcid a deplacer et renvoie une liste des deux AminoAcid a deplacer et dans quel
    direction si le deplacement est possible
    :param mat: matrice d'interet
    :param aa: AminoAcid d'interet
    :param prots: dictionnaire d'AminoAcid
    :return: 1 si deplacement impossible, une liste avec les deux acides amines a deplacer et la direction du
    deplacement si un deplacement est possible
    """
    aaname = aa.get_name()
    neighbors = check_neighbors(mat, aa, 1)
    neighbors_names = []
    neighl_names = []
    nlinkname = -2
    for neighbor in neighbors:
        neighbors_names.append(neighbor.get_name())
    neigh_linked = direct_neighbors(neighbors, aaname)
    neighborslink0 = check_neighbors(mat, neigh_linked[0], 1)
    neighborslink1 = check_neighbors(mat, neigh_linked[1], 1)
    if len(neighborslink0) == 0 or len(neighborslink1) == 0:
        return 1
    for nlink in neighborslink0:
        if nlink in neighborslink1 and nlink != aa:
            nlinkname = nlink.get_name()
            break
        else:
            return 1
    for neighl in neigh_linked:
        neighl_names.append(neighl.get_name())
    neighl_names.sort()
    if nlinkname == neighl_names[0] - 1:
        pred_limit = nlinkname
        succ_limit = neighl_names[1]
        movable_aa = nlinkname + 1
    elif nlinkname == neighl_names[1] + 1:
        pred_limit = neighl_names[0]
        succ_limit = nlinkname
        movable_aa = nlinkname - 1
    else:
        return 1
    coorpi, coorpj = prots[pred_limit].get_coor()
    coorsi, coorsj = prots[succ_limit].get_coor()
    if coorpi == coorsi:
        if mat[coorpi + 1][coorpj] == 0 and mat[coorsi + 1][coorsj] == 0:
            return [aa, prots[movable_aa], "down"]
        if mat[coorpi - 1][coorpj] == 0 and mat[coorsi - 1][coorsj] == 0:
            return [aa, prots[movable_aa], "up"]
    if coorpj == coorsj:
        if mat[coorpi][coorpj - 1] == 0 and mat[coorsi][coorsj - 1] == 0:
            return [aa, prots[movable_aa], "left"]
        if mat[coorsi][coorsj + 1] == 0 and mat[coorpi][coorpj + 1] == 0:
            return [aa, prots[movable_aa], "right"]
    return 1


def crank_move(mat, details, changei, changej):
    """
    Permet d'effectuer un mouvement dans une matrice selon les coordonnees a modifier
    :param mat: matrice d'interet
    :param details: une liste d'objet AminoAcid a deplacer dans la matrice
    :param changei: valeur de decalement de i dans la matrice
    :param changej: valeur de decalement de j dans la matrice
    :return:
    """
    aai, aaj = details[0].get_coor()
    maai, maaj = details[1].get_coor()
    mat[aai][aaj] = 0
    mat[maai][maaj] = 0
    mat[aai + changei][aaj + changej] = details[0]
    mat[maai + changei][maaj + changej] = details[1]
    details[0].update_coor(aai + changei, aaj + changej)
    details[1].update_coor(maai + changei, maaj + changej)


def crankshaft_move(mat, details):
    """
    Permet d'effectuer un deplacement de crankshaft en se basant sur les informations données par la fonction de test
    :param mat: matrice d'interet
    :param details: liste contenant les deux proteines a deplacer et la direction du deplacement
    :return: ne renvoie rien mais modifie directement la matrice et les coordonnees des AminoAcid
    """
    if details[2] == "down":
        crank_move(mat, details, 2, 0)
    if details[2] == "up":
        crank_move(mat, details, -2, 0)
    if details[2] == "left":
        crank_move(mat, details, 0, -2)
    if details[2] == "right":
        crank_move(mat, details, 0, 2)


def vshd_move(mat, aa, prots):
    """
    Fonction permettant d'effectuer un des different deplacement de type VSHD selon la possibilite a en effectuer un,
    ne fait rien si un deplacement VSHD n'est pas possible
    :param mat: matrice d'interet
    :param aa: AminoAcid d'interet
    :param prots: dictionnaire d'AminoAcid
    :return: renvoie un 0 si un mouvement aux extremites est effectue sinon ne renvoie rien et modifie directement la
    matrice si la modification est possible
    """
    neighbors = check_neighbors(mat, aa, 1)
    linked_n = direct_neighbors(neighbors, aa.get_name())
    if len(linked_n) == 1:
        vshd_move_ends(linked_n[0], mat, aa)
        return True
    kran = vshd_crankshaft_test(mat, aa, prots)
    if kran != 1:
        crankshaft_move(mat, kran)
        return True
    else:
        moved = vshd_move_base(neighbors, mat, aa)
        return moved


def energy_calc(mat, prots):
    """
    Fonction permettant de calculer l'energie de la structure etudie en se basant sur la matrice et les informations
    retrouve au sein des objets AminoAcid, ne sont etudie que les AminoAcid avec une polarite H (hydrophobe) en marquant
    les AminoAcid deja en interaction pour ne pas les reutiliser dans le calcul d'energie actuel
    :param mat: matrice d'interet
    :param prots: dictionnaire des AminoAcid
    :return: la valeur d'energie
    """
    hydrophobic = []
    energy = 0
    for key in prots:
        if prots[key].get_pol() == "H":
            hydrophobic.append(prots[key])
    for aah in hydrophobic:
        if aah.get_value() == 0:
            succ = aah.get_succ()
            pred = aah.get_pred()
            neighbors = check_neighbors(mat, aah, 1)
            for neighbor in neighbors:
                neiname = neighbor.get_name()
                neipol = neighbor.get_pol()
                neivalue = neighbor.get_value()
                if neiname != succ and neiname != pred and neipol == "H" and neivalue == 0:
                    energy -= 1
                    aah.set_value(1)
                    neighbor.set_value(1)
    return energy


def energy_reset(prots):
    """
    Fonction permettant de remettre le marqueur d'energie a 0 apres un calcul d'energie
    :param prots: dictionnaire des AminoAcid
    :return: rien, modifie directement les objets AminoAcid
    """
    for key in prots:
        prots[key].set_value(0)


def pull_move_prelim_test(aa, prots):
    """
    Permet de tester si trois AminoAcid sont alignés en prevision d'un pull-move en verifiant leurs coordonnees i et j
    respectifs
    :param aa: AminoAcid d'interet
    :param prots: dictionnaire d'AminoAcid
    :return: 0 si non aligne, 1 si aligne a l'horizontal et 2 si aligne a la vertical
    """
    coori, coorj = aa.get_coor()
    pred = aa.get_pred()
    succ = aa.get_succ()
    if succ == len(args.sequence[0]) or pred == 'none':
        return 0
    pcoori, pcoorj = prots[pred].get_coor()
    scoori, scoorj = prots[succ].get_coor()
    if coori == pcoori == scoori:
        return 1
    elif coorj == pcoorj == scoorj:
        return 2
    else:
        return 0


def pull_move_empty_test(mat, aa, side):
    """
    Permet de verifier si un pull-move est possible en verifiant s'il y a des espaces vides dans la matrice permettant
    ce deplacement et enregistre tous ces mouvement possible dans une liste.
    Un mouvement sera ensuite selectionnee de facon aleatoire dans cette liste est renvoye par la fonction si la liste
    n'est pas vide
    :param mat: matrice d'interet
    :param aa: AminoAcid d'interet
    :param side: 1 ou 2 selon le sens de deplacement disponible (resultat de la fonction pull_move_prelim_test())
    :return: 0 si pas de mouvement possible, coordonnees de deplacement ainsi que coordonnes du second AminoAcid a
    deplacer
    """
    coori, coorj = aa.get_coor()
    aaneighbors = check_neighbors(mat, aa, 0)
    down = mat[coori + 1][coorj]
    up = mat[coori - 1][coorj]
    left = mat[coori][coorj - 1]
    right = mat[coori][coorj + 1]
    potential_move = []
    if len(aaneighbors) == 0:
        return 0
    if side == 1:
        predd = mat[coori + 1][coorj - 1]
        predu = mat[coori - 1][coorj - 1]
        succd = mat[coori + 1][coorj + 1]
        succu = mat[coori - 1][coorj + 1]
        if predd == 0 and down == 0:
            potential_move.append([[coori + 1, coorj], [coori + 1, coorj - 1], [coori, coorj + 1]])
        if predu == 0 and up == 0:
            potential_move.append([[coori - 1, coorj], [coori - 1, coorj - 1], [coori, coorj + 1]])
        if succu == 0 and up == 0:
            potential_move.append([[coori - 1, coorj], [coori - 1, coorj + 1], [coori, coorj - 1]])
        if succd == 0 and down == 0:
            potential_move.append([[coori + 1, coorj], [coori + 1, coorj + 1], [coori, coorj - 1]])
    if side == 2:
        predl = mat[coori - 1][coorj - 1]
        predr = mat[coori - 1][coorj + 1]
        succl = mat[coori + 1][coorj - 1]
        succr = mat[coori + 1][coorj + 1]
        if predl == 0 and left == 0:
            potential_move.append([[coori, coorj - 1], [coori - 1, coorj - 1], [coori + 1, coorj]])
        if predr == 0 and right == 0:
            potential_move.append([[coori, coorj + 1], [coori - 1, coorj + 1], [coori + 1, coorj]])
        if succl == 0 and left == 0:
            potential_move.append([[coori, coorj - 1], [coori + 1, coorj - 1], [coori - 1, coorj]])
        if succr == 0 and right == 0:
            potential_move.append([[coori, coorj + 1], [coori + 1, coorj + 1], [coori - 1, coorj]])
    if len(potential_move) == 0:
        return 0
    else:
        move = random.choice(potential_move)
        return move


def pull_neighbor(i, j, mat):
    """
    Permet de trouver les voisins d'un acide amines en se basant sur les coordonnes et une matrice (specialement utilise
    dans les fonction de pull-move)
    :param i: coordonne i
    :param j: coordonne j
    :param mat: matrice d'interet
    :return: une liste de voisin (objet AminoAcid)
    """
    turn = [-1, 1]
    neigh = []
    for posturn in turn:
        if mat[i + posturn][j] != 0:
            neigh.append(mat[i + posturn][j])
        if mat[i][j + posturn] != 0:
            neigh.append(mat[i][j + posturn])
    return neigh


def pull_coor_update(temp, prots):
    """
    Permet de mettre a jour les coordonnees d'objets AminoAcid a la suite d'un pull_move reussi
    :param temp: dictionnaire de nom d'AminoAcid avec les nouvelles coordonnees
    :param prots: dictionnaire d'AminoAcid dont les coordonnees vont etre mis a jour
    :return: mets a jour
    """
    for aa in temp:
        prots[aa].update_coor(temp[aa][0], temp[aa][1])


def pull_move_try(path, neigh_test, current_coor, matbis, mat, previous_coor, prots, coor_temp, limit):
    """
    fonction permettant d'effectuer un pull-move jusqu'a sa finition (tirage) et de remplacer la matrice d'interet
    par la matrice miroir lorsque le mouvement est fini
    :param path: valeur permettant de determiner dans quelle sens le tirage se fait
    :param neigh_test: AminoAcid sur lequel le tirage se fait
    :param current_coor: coordonnees actuel
    :param matbis: matrice miroir ou les mouvement ont lieu jusqu'a la finition
    :param mat: matrice d'interet
    :param previous_coor: coordonnes precedente a conserver pour le prochain AminoAcid a tirer
    :param prots: dictionnaire d'AminoAcid
    :param coor_temp: dictionnaire ou les nouvelles coordonnees des AminoAcid sont conserves
    :param limit: nom d'AminoAcid limitant
    :return: 0 et 0 si le tirage est fini, current_coor et previous_coor si le tirage doit continuer
    """
    target = neigh_test[0]
    neighbors = pull_neighbor(current_coor[0], current_coor[1], matbis)
    if target == limit or (target + 1 in neighbors and target - 1 in neighbors):
        pull_coor_update(coor_temp, prots)
        mat = matbis
        del neigh_test[0]
        return 0, 0
    else:
        neigh_test.append(target + path)
        matbis[previous_coor[0]][previous_coor[1]] = prots[target + path]
        tcoori, tcoorj = prots[target + path].get_coor()
        matbis[tcoori][tcoorj] = 0
        coor_temp[target + path] = [previous_coor[0], previous_coor[1]]
        current_coor = [previous_coor[0], previous_coor[1]]
        ti, tj = prots[neigh_test[0]].get_coor()
        previous_coor = [ti, tj]
        del neigh_test[0]
        return current_coor, previous_coor


def pull_move(mat, aa, prots):
    """
    Fonction principal de tirage(pull-move), elle permet de tester si un mouvement de tirage est possible et si oui de
    l'effectuer et le mener a son terme
    :param mat: matrice d'interet
    :param aa: AminoAcid d'interet
    :param prots: dictionnaire d'AminoAcid
    :return: True si le mouvement a bien eu lieu et False si le mouvement n'est pas possible
    """
    coor_temp = {}
    matbis = mat
    preliminary_test = pull_move_prelim_test(aa, prots)
    if preliminary_test != 0:
        move_test = pull_move_empty_test(mat, aa, preliminary_test)
        if move_test != 0:
            secaa = matbis[move_test[2][0]][move_test[2][1]].get_name()
            aa_name = aa.get_name()
            coori, coorj = aa.get_coor()
            matbis[move_test[1][0]][move_test[1][1]] = aa
            matbis[coori][coorj] = 0
            matbis[move_test[0][0]][move_test[0][1]] = prots[secaa]
            scoori, scoorj = prots[secaa].get_coor()
            matbis[scoori][scoorj] = 0
            coor_temp[aa_name] = [move_test[1][0], move_test[1][1]]
            coor_temp[secaa] = [move_test[0][0], move_test[0][1]]
            previous_coor = [coori, coorj]
            current_coor = [coor_temp[secaa][0], coor_temp[secaa][1]]
            neigh_test = [secaa]
            while len(neigh_test) != 0:
                if secaa + 1 == aa_name:
                    current_coor, previous_coor = pull_move_try(-1, neigh_test, current_coor, matbis, mat,
                                                                previous_coor, prots, coor_temp,
                                                                prots[0].get_name())
                elif secaa - 1 == aa_name:
                    current_coor, previous_coor = pull_move_try(1, neigh_test, current_coor, matbis, mat, previous_coor,
                                                                prots, coor_temp,
                                                                len(args.sequence[0]) - 1)
            return True
    return False


def move(mat, aa, prots):
    """
    Fonction principale de deplacement qui permet de selectionner au hasard l'un des deux movement existant(vshd ou
    pull-move) et de l'effectuer et si celui selectionne n'est pas possible essaie de realiser l'autre mouvement
    :param mat: matrice d'interet
    :param aa: AminoAcid d'interet
    :param prots: dictionnaire d'AminoAcid d'interet
    :return: matrice d'interet modifie si mouvement possible sinon rien
    """
    move_choice = [0, 1]
    choice = random.choice(move_choice)
    if choice == 0:
        moved = pull_move(mat, aa, prots)
        if moved is False:
            vshd_move(mat, aa, prots)
    if choice == 1:
        moved = vshd_move(mat, aa, prots)
        if moved is False:
            pull_move(mat, aa, prots)


def retrieve_aas(prots):
    """
    Fonction permettant de recuperer les noms d'acide amine(leurs position) depuis un dictionnaire d'AminoAcid
    :param prots: dictionnaire d'AminoAcid
    :return: liste de nom/position d'AminoAcid
    """
    aas = []
    for key in prots:
        aas.append(key)
    return aas


def mcsearch(steps, mat, prots, aas, t=0.5):
    """
    Permet d'effectuer des repliement sur une sequence proteique en utilisant une méthode de Monte-Carlo base sur les
    interactions hydrophobes et les changements d'energie qu'elles entrainent.
    L'acide amine sur lequel une tentative de deplacement sera effectue est choisi de facon aleatoire
    :param steps: nombre d'etape de recherche
    :param mat: matrice d'interet
    :param prots: dictionnaire d'AminoAcid
    :param aas: liste de noms/positions d'AminoAcid
    :param t: temperature
    :return: liste des energies a chaque etape
    """
    energy_track = []
    for step in range(steps):
        matbis = mat
        aa_choice = random.choice(aas)
        aa = prots[aa_choice]
        move(matbis, aa, prots)
        energy = energy_calc(mat, prots)
        energy_reset(prots)
        energy_bis = energy_calc(matbis, prots)
        energy_reset(prots)
        energy_delta = energy_bis - energy
        if energy_delta <= 0:
            mat = matbis
            energy_track.append(energy)
        else:
            q = random.choice([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1])
            if q > math.exp(-energy_delta / t):
                mat = matbis
                energy_track.append(energy)
    return energy_track


def switch_aa(aa):
    """
    Switcher permettant de passer le code une lettre d'un acide amine a un code trois lettres
    :param aa: code une lettre d'un acide amine
    :return: code trois lettres d'un acide amine
    """
    switcher = {
        "D": "ASP",
        "E": "GLU",
        "A": "ALA",
        "R": "ARG",
        "N": "ASN",
        "C": "CYS",
        "Q": "GLN",
        "G": "GLY",
        "H": "HIS",
        "I": "ILE",
        "L": "LEU",
        "K": "LYS",
        "M": "MET",
        "F": "PHE",
        "P": "PRO",
        "S": "SER",
        "T": "THR",
        "W": "TRP",
        "Y": "TYR",
        "V": "VAL"
    }
    return switcher.get(aa, "UKN")


def pdb_printer(prots):
    """
    Permet de creer un fichier PDB du repliement effectuer a l'issue de MCsearch en utilisant les coordonnees et autres
    informations conserver au sein des objets AminoAcid
    :param prots: dictionnaire d'AminoAcid:
    :return: creer un document PDB
    """
    to_be_printed = []
    for aa in prots:
        at = "ATOM"
        na = "NA"
        x, y = prots[aa].get_coor()
        switched = switch_aa(prots[aa].get_hid_name())
        line = f"{at:6s}{prots[aa].get_name():5} {na:^4s} {switched:3s} A{0:4}    {x:8.3f}{y:8.3f}{0:8.3f}{0:6}{0:6}" \
               f"          {prots[aa].get_pol():>2s}\n"
        to_be_printed.append(line)
    with open("../results/myaafile.pdb", 'w') as filout:
        for line in to_be_printed:
            filout.write(line)


def plotter(track):
    """
    Permet d'afficher un graphique des variations d'energies en fonction de l'etape du MCsearch
    :param track: liste des energies a chaque etape
    :return: affiche un graphique
    """
    plt.plot(track)
    plt.xlabel("nombre d'etape du MCsearch")
    plt.ylabel("valeur d'energie")
    plt.title("Variation d'energie en fonction du nombre d'etape du MCsearch")
    plt.show()


def main():
    use = translation(args.sequence[0], AAapolaire)
    size = len(args.sequence[0])
    tab = mat_maker(size + 40)
    tab1 = mat_filler(tab, use)
    aa_names = retrieve_aas(use)
    printer(tab1)
    print("\n")
    if args.temperature:
        if args.step:
            step = int(args.step[0])
            energy_tracked = mcsearch(step, tab1, use, aa_names, args.temperature[0])
        else:
            energy_tracked = mcsearch(500, tab1, use, aa_names, args.temperature[0])
    else:
        if args.step:
            step = int(args.step[0])
            energy_tracked = mcsearch(step, tab1, use, aa_names)
        else:
            energy_tracked = mcsearch(500, tab1, use, aa_names)
    printer(tab1)
    pdb_printer(use)
    plotter(energy_tracked)


if __name__ == '__main__':
    main()
