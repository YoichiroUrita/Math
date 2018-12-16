Function SVD(o As Object)
'Sigular Value Decomposition
'refer from `Numerical recipes in C`
'https://www2.units.it/ipl/students_area/imm2/files/Numerical_Recipes.pdf
'6/12/2014

Dim rv1(), a(), w(), v() As Double
Dim g, scal, anorm, s, h, f As Double
Dim i, j, k, l, N, M As Integer

N = o.Columns.Count
M = o.Rows.Count

ReDim rv1(1 To N)
ReDim a(1 To M, 1 To N)
ReDim w(1 To N)
ReDim v(1 To N, 1 To N)

g = 0#
scal = 0#
anorm = 0#

For i = 1 To M
    For j = 1 To N
        a(i, j) = o(i, j)
    Next j
Next i

'Householder reduction to bidiagonal form
For i = 1 To N
    l = i + 1
    rv1(i) = scal * g
    g = 0#
    s = 0#
    scal = 0#
    
    If i <= M Then
        For k = i To M
            scal = scal + Abs(a(k, i))
        Next k
        If scal <> 0 Then
            For k = i To M
                a(k, i) = a(k, i) / scal
                s = s + a(k, i) * a(k, i)
            Next k
            
            f = a(i, i)
            If f <> 0 Then
                g = -(Sqr(s) * f / Abs(f))
            Else
                g = -Sqr(s)
            End If
            h = f * g - s
            a(i, i) = f - g
            For j = l To N
                s = 0#
                
                For k = i To M
                    s = s + a(k, i) * a(k, j)
                Next k
                f = s / h
                For k = i To M
                    a(k, j) = a(k, j) + f * a(k, i)
                Next k
                    
            Next j
            
            For k = i To M
                a(k, i) = a(k, i) * scal
            Next k
        End If
    End If
    
    w(i) = scal * g
    
    g = 0#
    s = 0#
    scal = 0#
    
    If i <= M And i <> N Then
        For k = l To N
            scal = scal + Abs(a(i, k))
        Next k
        If scal <> 0 Then
            For k = l To N
                a(i, k) = a(i, k) / scal
                s = s + a(i, k) * a(i, k)
            Next k
            f = a(i, l)
            If f = 0 Then
                g = -Sqr(s)
            Else
                g = -Sqr(s) * f / Abs(f)
            End If
            h = f * g - s
            a(i, l) = f - g
            For k = l To N
                rv1(k) = a(i, k) / h
            Next k
            For j = l To M
                s = 0#
                For k = l To N
                    s = s + a(j, k) * a(i, k)
                Next k
                For k = l To N
                    a(j, k) = a(j, k) + s * rv1(k)
                Next k
            Next j
            For k = l To N
                a(i, k) = a(i, k) * scal
            Next k
        End If
    End If
    If anorm < (Abs(w(i)) + Abs(rv1(i))) Then
        anorm = (Abs(w(i) + Abs(rv1(i))))
    End If
Next i

'accumulation of right-hand transformations
For i = N To 1 Step -1
    If i < N Then
        If g <> 0 Then
            'double division to avoid possible underflow
            For j = l To N
                v(j, i) = (a(i, j) / a(i, l)) / g
            Next j
            For j = l To N
                s = 0#
                For k = l To N
                    s = s + a(i, k) * v(k, j)
                Next k
                For k = l To N
                    v(k, j) = v(k, j) + s * v(k, i)
                Next k
            Next j
        End If
        For j = l To N
            v(i, j) = 0#
            v(j, i) = 0#
        Next j
    End If
    v(i, i) = 1#
    g = rv1(i)
    l = i
    If i = 1 Then
        Exit For
    End If
Next i

'accumulation of left-hand transformations
If M < N Then
    tmp = M
Else
    tmp = N
End If

For i = tmp To l Step -1
    l = i + 1
    g = w(i)
    For j = l To N
        a(i, j) = 0#
    Next j
    If g <> 0 Then
        g = 1# / g
        For j = l To N
            s = 0#
            For k = l To M
                s = s + a(k, i) * a(k, j)
            Next k
            f = (s / a(i, i)) * g
            For k = i To M
                a(k, j) = a(k, j) + f * a(k, i)
            Next k
        Next j
        For j = i To M
            a(j, i) = a(j, i) * g
        Next j
    Else
        For j = i To M
            a(j, i) = 0#
        Next j
    End If
    a(i, i) = a(i, i) + 1
    If i = l Then
        Exit For
    End If
Next i

'diagonalization of the bidiagonal form
':Loop over singular value, and over allowed iterations
For k = N To 1 Step -1
    For its = 1 To 30
        flag = 1
        For l = k To 1 Step -1
            nm = l - 1
            If Abs(rv1(l)) + anorm = anorm Then
                flag = 0
                Exit For
            End If
            If Abs(w(nm)) + anorm = anorm Then
                Exit For
            End If
            If l = 1 Then
                Exit For
            End If
        Next l
        'cancellation of rv1(l),if l>1
        If flag <> 0 Then
            c = 0#
            s = 0#
            For i = l To k
                f = s * rv1(i)
                rv1(i) = c * rv1(i)
                If Abs(f) + norm = norm Then
                    Exit For
                End If
                g = w(i)
                h = Sqr(f ^ 2 + g ^ 2)
                w(i) = h
                h = 1# / h
                c = g * h
                s = -f * g
                For j = 1 To M
                    y = a(j, nm)
                    z = a(j, i)
                    a(j, nm) = y * c + z * s
                    a(j, i) = z * c - y * s
                Next j
            Next i
        End If
        z = w(k)
        
        'convergence
        If l = k Then
            'singular value is made nonnegative
            If z < 0# Then
                w(k) = -z
                For j = 1 To N
                    v(j, k) = -v(j, k)
                Next j
            End If
            Exit For
        End If
        
        If its = 30 Then
            MsgBox ("no convergence in 30 svd iterations")
        End If
        
        'shift formbottom 2-by-2 minor
        x = w(l)
        nm = k - 1
        y = w(nm)
        g = rv1(nm)
        h = rv1(k)
        f = ((y - z) * (y + z) + (g + h)) / (2# * h * y)
        g = Sqr(f ^ 2 + 1#)
        f = ((x - z) * (x + z) + h * ((y / (f + Abs(g) * f / Abs(f))) - h)) / x
        'next qr transformation
        c = 1#
        s = 1#
        For j = 1 To nm
            i = j + 1
            g = rv1(i)
            y = w(i)
            h = s * g
            g = c * g
            z = Sqr(f ^ 2 + h ^ 2)
            rv1(j) = z
            c = f / z
            s = h / z
            f = x * c + g * s
            g = g * c - x * s
            h = y * s
            y = y * c
            For JJ = 1 To N
                x = v(JJ, j)
                z = v(JJ, i)
                v(JJ, j) = x * c + z * s
                v(JJ, i) = z * c - x * s
            Next JJ
            z = Sqr(f ^ 2 + h ^ 2)
            w(j) = z
            'rotation can be arbitrary if z=0
            If z <> 0 Then
                z = 1# / z
                c = f * z
                s = h * z
            End If
            f = c * g + s * y
            x = c * y - s * g
            For JJ = 1 To M
                y = a(JJ, j)
                z = a(JJ, i)
                a(JJ, j) = y * c + z * s
                a(JJ, i) = z * c - y * s
            Next JJ
        Next j
        rv1(l) = 0#
        rv1(k) = f
        w(k) = x
    Next its
    
    If k = 1 Then
        Exit For
    End If
Next k
    
SVD = w '[sigma]
' [M] = [U] * [sigma] * [V]^t
' where 
'[M] is m * n matrix 
'[U] is m * m matrix (left-singular vectors)
'[sigma] is diagonal m * n matrix (singular value)
'[V] is n * n matrix (right-singular vectors)
'
'if you want 
'[U], then use `a` instead of `w`
'[v]^t, then use `v` instead of `w`

End Function
