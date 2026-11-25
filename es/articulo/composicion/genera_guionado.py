import re, sys, glob, os
R, OUT = sys.argv[1:3]  # R: carpeta del artículo (con main.tex y secciones/); OUT: guionado.tex
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from silabas import hyphenated
EXCL = set('''abstract aligned axioms center checker choice document enumerate equation euler figure gather
hotmail itemize mailto mejiadilan multline pdfauthor pdfkeywords pdflang pdfsubject pdftitle print proof
propext queenon queenons section simkin sound subsec switchings table tabular zenodo lipschitz intro
doblecuenta preambulo articulo bibliografia codigo length minimum south style thick width comparacion
declaracion definicion discusion fisicos ilustracion interpretacion introduccion limite nucleo numerico
observacion proposicion reparacion verificacion exito vacios rcccccccccc rcccccccccccc
predictor constante insercion macro meso perfil preliminares resumen'''.split())
NAMES = set('''Nyström Nobel Agrawal Boyd Azuma Bowtell Cauchy Correia Dilan Euler Glock Gronwall Hoeffding Jensen Keevash Lipschitz Luria
Mathlib Mejia Munhá Nielsen Richardson Schwarz Simkin Stealth Steiner Steven Stirling Sudakov README AVISO
ArticuloA CierreDeLosSeisBloquesDelArticuloA DOIPerfil DOIcodigo DOIcodigoLean Classical
Bernoulli Bethe Bregman Eberhard Fukumizu Green Hilbert Ihara Jacobi Manners Mrazovi Neumann Newton Schmidt
Schur Slater Taylor Watanabe Vontobel Chertkov Chernyak Frobenius Perron Fredholm Carleman Weyl Ward Bass
Legendre Rodrigues Leibniz Kronecker Tonelli Bruns Poisson Laplace Szegő Kugelfunctionen Theorie Riemann
Gurvits Rademacher Stevens Schrijver Montanari Rizzo Yedidia Pitman'''.split())
files = [os.path.join(R, 'main.tex')] + sorted(glob.glob(R + '/secciones/*.tex')) + sorted(glob.glob(R + '/figuras/*.tex'))
files = [f for f in files if not f.endswith(('bibliografia.tex', 'simbolos.tex'))]
W = set()
for f in files:
    t = open(f, encoding='utf-8').read()
    t = re.sub(r'(?<!\\)%.*', '', t)
    t = re.sub(r'\\begin\{tikzpicture\}.*?\\end\{tikzpicture\}', ' ', t, flags=re.S)
    t = re.sub(r'\\(label|ref|eqref|cite|href|url|nolinkurl|input|texttt)\{[^}]*\}', ' ', t)
    t = re.sub(r'\\[A-Za-z]+\*?', ' ', t)
    for m in re.finditer(r'[A-Za-zÁÉÍÓÚÜÑáéíóúüñ]+', t):
        w = m.group(0)
        if len(w) < 5 or w in NAMES: continue
        w = w.lower()
        if w in EXCL: continue
        W.add(w)
items = sorted(h for h in (hyphenated(w) for w in W) if h and '-' in h)
ANCHO, MAXLIN = 108, 150   # líneas de a lo sumo 108 caracteres y archivos de a lo sumo 150 líneas
lines, cur = [], ''
for h in items:
    if len(cur) + len(h) + 1 > ANCHO:
        lines.append(cur); cur = h
    else:
        cur = (cur + ' ' + h).strip()
lines.append(cur)
hdr = ('% Puntos de corte de las palabras del artículo según la división silábica de la Ortografía de\n'
       '% la lengua española (RAE y ASALE, 2010); no se separan vocales contiguas. Sólo se usa cuando\n'
       '% la instalación de TeX no trae los patrones de guionado del español (véase preambulo.tex).\n')
cab = 5   # tres líneas de comentario, la apertura y el cierre
if len(lines) + cab <= MAXLIN:
    open(OUT, 'w', encoding='utf-8').write(hdr + '\\babelhyphenation[spanish]{\n' + '\n'.join(lines) + '\n}\n')
else:
    # Se reparte en dos archivos; el primero carga el segundo.
    n1 = MAXLIN - cab - 1
    OUTB = OUT[:-4] + '_b.tex'
    open(OUT, 'w', encoding='utf-8').write(hdr + '\\babelhyphenation[spanish]{\n' + '\n'.join(lines[:n1])
                                         + '\n}\n\\input{' + os.path.basename(OUTB)[:-4] + '}\n')
    open(OUTB, 'w', encoding='utf-8').write('% Continuación de guionado.tex.\n\\babelhyphenation[spanish]{\n'
                                          + '\n'.join(lines[n1:]) + '\n}\n')
print(len(W), 'palabras;', len(items), 'con cortes;', len(lines), 'lineas')
