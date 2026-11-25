"""Division silabica del espanol para puntos de guionado (reglas de la Ortografia de la RAE, 2010),
en version conservadora: nunca se separan vocales contiguas ni digrafos (ch, ll, rr)."""
import re, sys
VOW = set('aeiouáéíóúü')
ONSET = {'pl','pr','bl','br','fl','fr','cl','cr','gl','gr','dr','tr','kl','kr'}

def units(w):
    out, i = [], 0
    while i < len(w):
        two = w[i:i+2]
        if two in ('ch', 'll', 'rr'):
            out.append(('C', two)); i += 2; continue
        c = w[i]
        if c == 'y':
            nxt = w[i+1] if i+1 < len(w) else ''
            out.append(('C' if nxt in VOW else 'V', c)); i += 1; continue
        out.append(('V' if c in VOW else 'C', c)); i += 1
    return out

def breaks(word):
    """Posiciones (indices de caracter) donde se permite cortar la palabra en minusculas."""
    w = word.lower()
    if not re.fullmatch(r'[a-zñáéíóúü]+', w) or 'w' in w:
        return None
    if re.match(r'sub[lr]', w):          # sub+l/r: el corte silabico no coincide con el prefijo
        return None
    u = units(w)
    pos, starts = 0, []
    for t, s in u:
        starts.append(pos); pos += len(s)
    vidx = [k for k, (t, s) in enumerate(u) if t == 'V']
    res = []
    for a, b in zip(vidx, vidx[1:]):
        cons = list(range(a+1, b))
        if not cons:
            continue                      # vocales contiguas: no se separan
        if len(cons) == 1:
            cut = cons[0]
        else:
            p, q = u[cons[-2]][1], u[cons[-1]][1]
            if p + q == 'tl':
                continue                  # tl: division regional variable, se evita
            cut = cons[-2] if (p + q) in ONSET else cons[-1]
        res.append(starts[cut])
    n = len(w)
    return [k for k in res if k >= 2 and n - k >= 2]

def hyphenated(word):
    b = breaks(word)
    if not b:
        return None
    out, last = [], 0
    for k in b:
        out.append(word[last:k]); last = k
    out.append(word[last:])
    return '-'.join(out)

if __name__ == '__main__':
    for w in sys.argv[1:]:
        print(w, '->', hyphenated(w))
