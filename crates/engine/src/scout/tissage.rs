//! Clayonnage et PH (géométrie pure, sans le noyau) : les claies tressées entre des supports et les
//! cordes tendues entre deux traverses, comme dans CathIA. Le résultat est une liste de tronçons
//! droits (un cylindre chacun) ; la commande les construit et les range dans un composant.

use super::pole::{Build, JEU, Op, Pole, V3, add, cross, dot, len, mul, norm, sub};

/// Un tronçon droit de claie ou de corde, de `a` à `b`, de diamètre `d`.
#[derive(Clone, Debug, PartialEq)]
pub struct Troncon {
    pub a: V3,
    pub b: V3,
    pub d: f64,
}

/// Une claie : la spline qui passe par `points` (un bout, chaque support, l'autre bout), de diamètre `d`.
#[derive(Clone, Debug, PartialEq)]
pub struct Claie {
    pub points: Vec<V3>,
    pub d: f64,
}

/// Les réglages du clayonnage (millimètres), ceux de CathIA.
#[derive(Clone, Debug)]
pub struct Clayonnage {
    /// Diamètre des claies (8 à 20 mm).
    pub d: f64,
    /// Pas entre deux claies ; 0 = jointives (pas = diamètre).
    pub pas: f64,
    /// Marge laissée à chaque extrémité des supports.
    pub marge: f64,
    /// Début du clayage après la marge.
    pub debut: f64,
    /// Longueur clayée ; 0 = jusqu'au bout.
    pub zone: f64,
    /// Débord au-delà des supports extrêmes (5 mm au plus).
    pub debord: f64,
    /// La première claie passe sur le premier support (sinon dessous).
    pub dessus: bool,
}

impl Default for Clayonnage {
    fn default() -> Self {
        Clayonnage { d: 20.0, pas: 0.0, marge: 10.0, debut: 0.0, zone: 0.0, debord: 0.0, dessus: false }
    }
}

/// Au plus autant de claies (au-delà, la construction serait trop longue : augmenter le pas).
pub const MAX_CLAIES: usize = 300;

/// Le point de l'axe de `p` le plus proche de `q`.
fn sur_axe(p: &Pole, q: V3) -> V3 {
    add(p.origin, mul(p.dir, dot(sub(q, p.origin), p.dir)))
}

/// Les claies tressées entre des supports parallèles (au moins deux, rangés dans l'ordre où les claies
/// les traversent). Chaque claie passe alternativement dessus et dessous ; une claie sur deux commence
/// de l'autre côté, si bien que deux claies voisines se croisent sur chaque support.
pub fn clayonnage(supports: &[Pole], o: &Clayonnage) -> Result<Vec<Claie>, String> {
    if supports.len() < 2 {
        return Err("choisissez au moins deux supports (des traverses parallèles)".into());
    }
    if !(8.0..=40.0).contains(&o.d) {
        return Err("le diamètre des claies va de 8 à 40 mm".into());
    }
    let u = supports[0].dir;
    if supports.iter().any(|p| dot(p.dir, u).abs() < 0.985) {
        return Err("les supports doivent être parallèles".into());
    }
    // la direction des claies : du premier support vers le dernier, perpendiculaire aux supports
    let o0 = supports[0].origin;
    let mut ordre: Vec<(f64, &Pole)> = Vec::new();
    let ecart = |p: &Pole| {
        let q = sur_axe(p, o0);
        sub(q, o0)
    };
    let mut v = V3::default();
    for p in supports.iter().skip(1) {
        let e = ecart(p);
        if len(e) > len(v) {
            v = e;
        }
    }
    if len(v) < o.d {
        return Err("les supports sont confondus : écartez-les".into());
    }
    let v = norm(v);
    for p in supports {
        ordre.push((dot(ecart(p), v), p));
    }
    ordre.sort_by(|a, b| a.0.total_cmp(&b.0));
    if ordre.windows(2).any(|w| w[1].0 - w[0].0 < o.d) {
        return Err("deux supports sont trop proches l'un de l'autre".into());
    }
    // le côté « dessus » : perpendiculaire aux supports et aux claies, vers le haut
    let mut w = norm(cross(u, v));
    if w[2] < 0.0 {
        w = mul(w, -1.0);
    }
    // la partie commune des supports, le long de u (mesurée depuis o0)
    let (mut t0, mut t1) = (f64::NEG_INFINITY, f64::INFINITY);
    for (_, p) in &ordre {
        let (a, b) = (dot(sub(p.p_start(), o0), u), dot(sub(p.end(), o0), u));
        t0 = t0.max(a.min(b));
        t1 = t1.min(a.max(b));
    }
    let debut = t0 + o.marge.max(0.0) + o.debut.max(0.0);
    let fin = if o.zone > 0.0 { (debut + o.zone).min(t1 - o.marge.max(0.0)) } else { t1 - o.marge.max(0.0) };
    if fin - debut < o.d {
        return Err("les supports ne se recouvrent pas assez pour poser une claie".into());
    }
    let pas = if o.pas > 0.0 { o.pas.max(o.d) } else { o.d };
    let n = ((fin - debut) / pas).floor() as usize + 1;
    if n > MAX_CLAIES {
        return Err(format!("{n} claies : plus que les {MAX_CLAIES} permises (augmentez le pas ou réduisez la longueur clayée)"));
    }
    let debord = o.debord.clamp(0.0, 5.0);
    let mut claies = Vec::with_capacity(n);
    for i in 0..n {
        let t = debut + pas * i as f64 + o.d / 2.0;
        if t > fin + 1e-6 {
            break;
        }
        let repere = add(o0, mul(u, t));
        let points: Vec<V3> = ordre
            .iter()
            .enumerate()
            .map(|(k, (_, p))| {
                let dessus = ((i + k) % 2 == 0) == o.dessus;
                let c = sur_axe(p, repere);
                add(c, mul(w, if dessus { p.radius() + o.d / 2.0 } else { -(p.radius() + o.d / 2.0) }))
            })
            .collect();
        // la claie dépasse un peu les supports extrêmes (d/2 pour couvrir le support, plus le débord),
        // à plat : la spline qui passe par ces points ondule sur et sous les supports
        let dep = o.d / 2.0 + debord;
        let mut pts = Vec::with_capacity(points.len() + 2);
        if let (Some(first), Some(last)) = (points.first(), points.last()) {
            pts.push(add(*first, mul(v, -dep)));
            pts.extend(points.iter().copied());
            pts.push(add(*last, mul(v, dep)));
        }
        claies.push(Claie { points: pts, d: o.d });
    }
    Ok(claies)
}

/// Le PH (paraboloïde hyperbolique, décoratif) : `n` cordes de diamètre `d` entre deux traverses
/// qui se rejoignent. Le point de départ est la jonction des deux axes ; chaque corde relie un point
/// de la traverse 1 à un point de la traverse 2, l'un s'éloignant de la jonction pendant que l'autre
/// s'en rapproche. `la`, `lb` : longueurs utilisées sur chaque traverse (0 = jusqu'au bout).
pub fn ph(a: &Pole, b: &Pole, n: usize, la: f64, lb: f64, d: f64) -> Result<Vec<Troncon>, String> {
    if !(3..=80).contains(&n) {
        return Err("de 3 à 80 cordes".into());
    }
    if !(1.0..=20.0).contains(&d) {
        return Err("le diamètre de la corde va de 1 à 20 mm".into());
    }
    let (u, v) = (a.dir, b.dir);
    let uv = dot(u, v);
    if uv.abs() > 0.995 {
        return Err("les deux traverses sont parallèles : un PH demande deux traverses qui se croisent ou se rejoignent".into());
    }
    // points les plus proches entre les deux axes
    let w0 = sub(a.origin, b.origin);
    let (d1, e1) = (dot(u, w0), dot(v, w0));
    let den = 1.0 - uv * uv;
    let (sa, sb) = ((uv * e1 - d1) / den, (e1 - uv * d1) / den);
    let (ja, jb) = (add(a.origin, mul(u, sa)), add(b.origin, mul(v, sb)));
    if len(sub(ja, jb)) > a.radius() + b.radius() + 50.0 {
        return Err("les deux traverses ne se rejoignent pas (écart de plus de 5 cm entre elles)".into());
    }
    // sur chaque traverse : vers le bout le plus éloigné de la jonction
    let cote = |p: &Pole, j: V3| {
        let (d0, d1) = (len(sub(p.p_start(), j)), len(sub(p.end(), j)));
        if d0 > d1 { (norm(sub(p.p_start(), j)), d0) } else { (norm(sub(p.end(), j)), d1) }
    };
    let ((ua, maxa), (ub, maxb)) = (cote(a, ja), cote(b, jb));
    let la = if la > 0.0 { la.min(maxa) } else { maxa };
    let lb = if lb > 0.0 { lb.min(maxb) } else { maxb };
    if la < 100.0 || lb < 100.0 {
        return Err("les traverses sont trop courtes au-delà de leur jonction (10 cm au moins)".into());
    }
    let mut cordes = Vec::with_capacity(n);
    for i in 1..=n {
        let f = i as f64 / (n + 1) as f64;
        let (p, q) = (add(ja, mul(ua, la * f)), add(jb, mul(ub, lb * (1.0 - f))));
        let dir = sub(q, p);
        if len(dir) < a.radius() + b.radius() + 10.0 {
            continue;
        }
        let dir = norm(dir);
        // la corde part de la surface des traverses, pas de leur axe
        cordes.push(Troncon { a: add(p, mul(dir, a.radius())), b: sub(q, mul(dir, b.radius())), d });
    }
    Ok(cordes)
}

/// Les plateaux posés dans des quarts de bois (CathIA) : claies droites ou planches de rondins fendus.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Garniture {
    Claies,
    Planches,
}

/// Les réglages d'un plateau posé dans des quarts de bois (millimètres).
#[derive(Clone, Debug)]
pub struct Plateau {
    pub garniture: Garniture,
    /// Diamètre des claies.
    pub d: f64,
    /// Épaisseur et largeur des planches.
    pub ep: f64,
    pub largeur: f64,
    /// Jeu entre deux éléments.
    pub jeu: f64,
    /// Début après l'extrémité des supports, et fin depuis l'autre extrémité (0 = jusqu'au bout).
    pub debut: f64,
    pub fin: f64,
    /// Tailler les quarts de bois dans les supports (sinon, éléments posés dessus).
    pub quart: bool,
    /// Cheviller chaque planche dans chaque support, et le diamètre des chevilles.
    pub chevilles: bool,
    pub dp: f64,
}

impl Default for Plateau {
    fn default() -> Self {
        Plateau {
            garniture: Garniture::Claies,
            d: 20.0,
            ep: 25.0,
            largeur: 100.0,
            jeu: 0.0,
            debut: 0.0,
            fin: 0.0,
            quart: true,
            chevilles: true,
            dp: 10.0,
        }
    }
}

/// Au plus autant d'éléments sur un plateau.
pub const MAX_ELEMENTS: usize = 200;

/// Pose un plateau entre les supports `a` et `b` d'un assemblage en cours : un quart de bois taillé dans
/// chaque support (le quartier du haut, côté intérieur), puis les claies ou les planches posées en
/// travers dans ces quarts de bois, chaque planche chevillée dans chaque support si demandé.
/// Les supports doivent être parallèles et horizontaux. Renvoie le nombre d'éléments posés.
pub fn plateau(b: &mut Build, ia: usize, ib: usize, o: &Plateau) -> Result<usize, String> {
    let (pa, pb) = match (b.poles.get(ia), b.poles.get(ib)) {
        (Some(x), Some(y)) if ia != ib => (x.clone(), y.clone()),
        _ => return Err("choisissez deux supports différents".into()),
    };
    let nombres = [o.d, o.ep, o.largeur, o.jeu, o.debut, o.fin, o.dp];
    if nombres.iter().any(|v| !v.is_finite() || *v < 0.0) {
        return Err("les dimensions sont des nombres positifs".into());
    }
    let u = pa.dir;
    if dot(pb.dir, u).abs() < 0.985 {
        return Err("les deux supports doivent être parallèles".into());
    }
    if u[2].abs() > 0.1 || pb.dir[2].abs() > 0.1 {
        return Err("un plateau dans des quarts de bois se pose sur des supports horizontaux".into());
    }
    let ecart = sub(sur_axe(&pb, pa.origin), pa.origin);
    if len(ecart) < pa.radius() + pb.radius() + 20.0 {
        return Err("les deux supports sont trop proches".into());
    }
    let v = norm([ecart[0], ecart[1], 0.0]);
    let haut: V3 = [0.0, 0.0, 1.0];
    // la partie commune des supports, le long de u
    let (a0, a1) = (dot(sub(pa.p_start(), pa.origin), u), dot(sub(pa.end(), pa.origin), u));
    let (b0, b1) = (dot(sub(pb.p_start(), pa.origin), u), dot(sub(pb.end(), pa.origin), u));
    let (t0, t1) = (a0.min(a1).max(b0.min(b1)), a0.max(a1).min(b0.max(b1)));
    let debut = t0 + o.debut;
    let fin = if o.fin > 0.0 { t1 - o.fin } else { t1 };
    let zone = fin - debut;
    let (pas, epaisseur) = match o.garniture {
        Garniture::Claies => {
            if !(8.0..=40.0).contains(&o.d) {
                return Err("le diamètre des claies va de 8 à 40 mm".into());
            }
            (o.d + o.jeu, o.d)
        }
        Garniture::Planches => {
            if !(10.0..=80.0).contains(&o.ep) || !(40.0..=400.0).contains(&o.largeur) {
                return Err("planches : épaisseur de 10 à 80 mm, largeur de 40 à 400 mm".into());
            }
            (o.largeur + o.jeu, o.ep)
        }
    };
    let largeur = if o.garniture == Garniture::Claies { o.d } else { o.largeur };
    if zone < largeur {
        return Err("la zone à garnir est plus courte qu'un élément".into());
    }
    let n = ((zone - largeur) / pas).floor() as usize + 1;
    if n > MAX_ELEMENTS {
        return Err(format!("{n} éléments : plus que les {MAX_ELEMENTS} permis (augmentez le jeu)"));
    }
    // le quart de bois : le quartier du haut, côté intérieur, sur toute la zone
    let quart = |p: &Pole, interieur: V3| -> Result<Op, String> {
        let (xl, yl) = (p.local_dir(haut), p.local_dir(interieur));
        if xl[0].abs() < 0.95 {
            return Err(format!("« {} » est tourné sur lui-même : le quart de bois demande une pièce posée à plat", p.role));
        }
        let r = p.radius();
        let milieu = add(sur_axe(p, add(pa.origin, mul(u, (debut + fin) / 2.0))), add(mul(haut, r / 2.0 + 1.0), mul(interieur, r / 2.0 + 1.0)));
        let c = p.local(milieu);
        let _ = yl;
        Ok(Op::Notch { c, h: [r / 2.0 + 1.0, r / 2.0 + 1.0, zone / 2.0 + JEU] })
    };
    // la hauteur où reposent les éléments : l'axe des supports dans un quart de bois, leur dessus sinon
    let (za, zb) = if o.quart { (0.0, 0.0) } else { (pa.radius(), pb.radius()) };
    if o.quart {
        let qa = quart(&pa, v)?;
        let qb = quart(&pb, mul(v, -1.0))?;
        b.poles[ia].push(qa)?;
        b.poles[ib].push(qb)?;
    }
    for i in 0..n {
        let t = debut + largeur / 2.0 + pas * i as f64;
        let repere = add(pa.origin, mul(u, t));
        let (ca, cb) = (sur_axe(&pa, repere), sur_axe(&pb, repere));
        let (p, q) = (add(ca, mul(haut, za + epaisseur / 2.0)), add(cb, mul(haut, zb + epaisseur / 2.0)));
        let element = match o.garniture {
            Garniture::Claies => Pole::bar("claie", o.d, p, q),
            Garniture::Planches => Pole::board("planche", o.ep, o.largeur, p, q, haut),
        };
        b.poles.push(element);
        let k = b.poles.len() - 1;
        if o.garniture == Garniture::Planches && o.chevilles {
            if !(4.0..=30.0).contains(&o.dp) {
                return Err("le diamètre des chevilles va de 4 à 30 mm".into());
            }
            // une cheville verticale par support, au milieu de la partie posée sur lui
            for (support, centre, r) in [(ia, ca, pa.radius()), (ib, cb, pb.radius())] {
                let vers = norm(sub(if support == ia { cb } else { ca }, centre));
                let pied = add(centre, mul(vers, r / 2.0));
                let total = o.ep + r;
                let haut_cheville = add(pied, mul(haut, (if o.quart { 0.0 } else { r }) + o.ep));
                let k_dir = mul(haut, -1.0);
                b.peg_line(k, support, sub(haut_cheville, mul(haut, total / 2.0)), k_dir, o.dp, total)?;
            }
        }
    }
    Ok(n)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn traverse(y: f64, z: f64) -> Pole {
        Pole::bar("support", 60.0, [0.0, y, z], [1000.0, y, z])
    }

    #[test]
    fn les_claies_alternent_dessus_et_dessous() {
        let s = [traverse(0.0, 500.0), traverse(300.0, 500.0), traverse(600.0, 500.0)];
        let c = clayonnage(&s, &Clayonnage::default()).expect("clayonnage");
        assert!(c.len() > 40 && c.len() < 51, "{} claies jointives de 20 mm sur 980 mm", c.len());
        for claie in &c {
            assert_eq!(claie.points.len(), 5, "un bout, trois supports, l'autre bout");
        }
        // première claie : sous le premier support, dessus le deuxième, sous le troisième
        let z: Vec<f64> = c[0].points[1..4].iter().map(|p| p[2]).collect();
        assert!(z[0] < 500.0 && z[1] > 500.0 && z[2] < 500.0, "{z:?}");
        // la claie voisine fait l'inverse
        assert!(c[1].points[1][2] > 500.0 && c[1].points[2][2] < 500.0);
        // les claies touchent les supports : centre à rayon du support + rayon de la claie
        assert!((c[0].points[2][2] - (500.0 + 30.0 + 10.0)).abs() < 1e-6);
        // les bouts dépassent les supports extrêmes à plat
        assert!((c[0].points[0][2] - c[0].points[1][2]).abs() < 1e-9 && c[0].points[0][1] < 0.0);
    }

    #[test]
    fn l_ordre_des_supports_ne_compte_pas() {
        let a = [traverse(0.0, 500.0), traverse(600.0, 500.0), traverse(300.0, 500.0)];
        let b = [traverse(0.0, 500.0), traverse(300.0, 500.0), traverse(600.0, 500.0)];
        assert_eq!(clayonnage(&a, &Clayonnage::default()), clayonnage(&b, &Clayonnage::default()));
    }

    #[test]
    fn le_pas_la_marge_et_la_zone_sont_respectes() {
        let s = [traverse(0.0, 500.0), traverse(400.0, 500.0)];
        let o = Clayonnage { pas: 50.0, marge: 100.0, zone: 300.0, ..Clayonnage::default() };
        let c = clayonnage(&s, &o).expect("clayonnage");
        assert_eq!(c.len(), 6, "300 mm au pas de 50 : 6 claies de 20 mm entières");
        assert!((c[0].points[0][0] - 110.0).abs() < 1e-6, "première claie à marge + rayon");
    }

    #[test]
    fn les_cas_impossibles_sont_refuses() {
        let o = Clayonnage::default();
        assert!(clayonnage(&[traverse(0.0, 500.0)], &o).is_err(), "un seul support");
        let croise = Pole::bar("support", 60.0, [500.0, -200.0, 500.0], [500.0, 800.0, 500.0]);
        assert!(clayonnage(&[traverse(0.0, 500.0), croise], &o).is_err(), "supports non parallèles");
        assert!(clayonnage(&[traverse(0.0, 500.0), traverse(0.0, 500.0)], &o).is_err(), "confondus");
        let fin = Clayonnage { pas: 1.0, d: 8.0, ..Clayonnage::default() };
        let long = [Pole::bar("s", 60.0, [0.0, 0.0, 0.0], [5000.0, 0.0, 0.0]), Pole::bar("s", 60.0, [0.0, 300.0, 0.0], [5000.0, 300.0, 0.0])];
        assert!(clayonnage(&long, &fin).is_err(), "trop de claies");
        assert!(clayonnage(&long, &Clayonnage { d: f64::NAN, ..Clayonnage::default() }).is_err(), "NaN");
    }

    #[test]
    fn le_ph_tend_ses_cordes_entre_les_deux_traverses() {
        let a = Pole::bar("a", 60.0, [0.0, 0.0, 0.0], [1000.0, 0.0, 0.0]);
        let b = Pole::bar("b", 60.0, [0.0, 0.0, 60.0], [0.0, 0.0, 1060.0]);
        let c = ph(&a, &b, 10, 0.0, 0.0, 4.0).expect("ph");
        assert_eq!(c.len(), 10);
        // la première corde part près de la jonction sur a et loin sur b
        assert!(c[0].a[0] < 200.0 && c[0].b[2] > 800.0);
        assert!(c[9].a[0] > 800.0 && c[9].b[2] < 300.0);
        assert!(ph(&a, &Pole::bar("p", 60.0, [0.0, 500.0, 0.0], [1000.0, 500.0, 0.0]), 10, 0.0, 0.0, 4.0).is_err(), "parallèles");
        assert!(ph(&a, &Pole::bar("p", 60.0, [0.0, 900.0, 0.0], [0.0, 900.0, 1000.0]), 10, 0.0, 0.0, 4.0).is_err(), "trop loin");
        assert!(ph(&a, &b, 1, 0.0, 0.0, 4.0).is_err(), "trop peu de cordes");
    }

    #[test]
    fn un_plateau_de_claies_se_pose_dans_des_quarts_de_bois() {
        let mut b = Build { poles: vec![traverse(0.0, 500.0), traverse(600.0, 500.0)], ..Default::default() };
        let n = plateau(&mut b, 0, 1, &Plateau::default()).expect("plateau");
        assert_eq!(n, 50, "1000 mm de claies jointives de 20 mm");
        assert!(matches!(b.poles[0].ops.last(), Some(Op::Notch { .. })) && matches!(b.poles[1].ops.last(), Some(Op::Notch { .. })));
        // les claies reposent sur l'axe des supports (fond du quart de bois)
        assert!((b.poles[2].origin[2] - 510.0).abs() < 1e-6);
        assert!((len(sub(b.poles[2].end(), b.poles[2].p_start())) - 600.0).abs() < 1e-6);
    }

    #[test]
    fn des_planches_chevillees_dans_chaque_support() {
        let mut b = Build { poles: vec![traverse(0.0, 500.0), traverse(600.0, 500.0)], ..Default::default() };
        let o = Plateau { garniture: Garniture::Planches, jeu: 3.0, ..Plateau::default() };
        let n = plateau(&mut b, 0, 1, &o).expect("planches");
        assert_eq!(n, 9, "1000 mm de planches de 100 mm avec 3 mm de jeu");
        assert_eq!(b.poles.len(), 2 + n * 3, "chaque planche et ses deux chevilles");
        assert!(b.poles[0].ops.iter().filter(|o| matches!(o, Op::Hole { .. })).count() == n, "un trou par planche dans le support");
        let mut c =
            Build { poles: vec![traverse(0.0, 500.0), Pole::bar("p", 60.0, [0.0, 300.0, 500.0], [0.0, 300.0, 1500.0])], ..Default::default() };
        assert!(plateau(&mut c, 0, 1, &o).is_err(), "supports non parallèles");
        let mut d = Build { poles: vec![traverse(0.0, 500.0), traverse(600.0, 500.0)], ..Default::default() };
        assert!(plateau(&mut d, 0, 1, &Plateau { d: f64::NAN, ..Plateau::default() }).is_err());
    }
}
