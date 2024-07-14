      ! VASP subroutine
      FUNCTION LENGTH(STRING)
      USE mod_comp
      IMPLICIT REAL(DP) (A-H,O-Z)
! Returns the position of the last non-blank character in STRING
      CHARACTER*(*)   STRING
      CHARACTER*256   B8
      CHARACTER*128   B7
      CHARACTER*64    B6
      CHARACTER*32    B5
      INTEGER         LENGTH,LEN,L,I

      SAVE            B5,B6,B7,B8
      DATA            B5 /' '/, B6 /' '/, B7 /' '/, B8 /' '/


      L=LEN(STRING)
! Very crude 'scan' for 'order of length' (performance!!) ...:
    6 IF (L.GE.256) THEN
         IF (STRING(L-255:L).EQ.B8) THEN
            L=L-256
            GOTO 7
         ELSE IF (STRING(L-127:L).EQ.B7) THEN
            L=L-128
            GOTO 7
         ELSE IF (STRING(L-63:L).EQ.B6) THEN
            L=L-64
            GOTO 7
         ELSE IF (STRING(L-31:L).EQ.B5) THEN
            L=L-32
            GOTO 7
         ENDIF
      ENDIF
    7 IF (L.GE.128) THEN
         IF (STRING(L-127:L).EQ.B7) THEN
            L=L-128
            GOTO 8
         ELSE IF (STRING(L-63:L).EQ.B6) THEN
            L=L-64
            GOTO 8
         ELSE IF (STRING(L-31:L).EQ.B5) THEN
            L=L-32
            GOTO 8
         ENDIF
      ENDIF
    8 IF (L.GE.64) THEN
         IF (STRING(L-63:L).EQ.B6) THEN
            L=L-64
            GOTO 9
         ELSE IF (STRING(L-31:L).EQ.B5) THEN
            L=L-32
            GOTO 9
         ENDIF
      ENDIF
    9 IF (L.GE.32) THEN
         IF (STRING(L-31:L).EQ.B5) L=L-32
      ENDIF
! Here we should have reached some point where either L is much
! smaller than LEN(STRING)/very small at all or more general where
! L is quite close to the result for length (so that all goes very
! quick now ... --- function LENGTH should show high performance).
      LENGTH=0
      DO 10 I=L,1,-1
         IF (STRING(I:I).NE.' ') THEN
            LENGTH=I
            GOTO 11
         ENDIF
   10 CONTINUE
   11 RETURN
      END


      FUNCTION NWORDS(STRING)
      USE mod_comp
      IMPLICIT REAL(DP) (A-H,O-Z)
! Find the number of blank-delimited words within some string:
      CHARACTER*(*) STRING
      LOGICAL BLANK
      INTEGER NWORDS,L,LENGTH,I
      EXTERNAL LENGTH

      L=LENGTH(STRING)
      NWORDS=0
      BLANK=.TRUE.
      DO 10 I=1,L
         IF (BLANK.AND.(STRING(I:I).NE.' ')) THEN
            BLANK=.FALSE.
            NWORDS=NWORDS+1
         ELSE IF ((.NOT.BLANK).AND.(STRING(I:I).EQ.' ')) THEN
            BLANK=.TRUE.
         ENDIF
   10 CONTINUE
      RETURN
      END


      FUNCTION NITEMS(STRING,WORK,EXPAND,TYPE)
      USE mod_comp
      IMPLICIT REAL(DP) (A-H,O-Z)
! Extract how many data items had to be read if one wanted to read
! numerical data from a STRING ... . --> It is not simply the number
! of words (because such constructs like 120*0.2 -- counting as 120
! different data in this example -- are also allowed in FORTRAN!!).
! If one wishes one can resolve such constructs using EXPAND=.TRUE.
! (i.e. STRING will be changed so that no more such constructs occur
! by translating it in a long list of single items -- of course here
! some restrictions apply: STRING and WORK must be long enough to
! hold the result and the single items may not exceed 255 characters
! and it should be noted that this acts also implicitly like STRIP,
! i.e. leading blanks are deleted and all items are separated by a
! single blank after all ... . And additionally one can also specify
! 'type verification' (means entering some given data via TYPE being
!   - 'L'   logical items only
!   - 'I'   integer numbers only
!   - 'F'   floating point items only
!   - 'C'   complex numbers only
!   - or any other type suppressing the check ...
! will result in a validity-check of the items found, breaking the
! searching and counting at the first invalid item found ...).
      CHARACTER*(*) STRING,WORK
      CHARACTER*1   CHECK,TYPE



      CHARACTER*255 BUFFER,FORM,NUMBER,DUMMY

      LOGICAL EXPAND
      INTEGER NITEMS,NWORDS,NOCCUR,LS,LW,L0,L,I,N,IC,IMULT,LENGTH,J,LEN
      INTRINSIC LEN
      EXTERNAL NWORDS,NOCCUR,LENGTH

      LS=LEN(STRING)
      LW=LEN(WORK)
      IF (EXPAND) CALL STRIP(STRING,L,'B')
  100 NITEMS=0
      N=NWORDS(STRING)
      DO 400 I=1,N
! scan word by word ...
         CALL SUBWRD(STRING,WORK,I,1)
         IC=NOCCUR(WORK,'*',0)
! invalid item here ---> stop counting items and say good bye ... !
         IF (IC.GT.1) GOTO 500
         IF (IC.EQ.0) THEN
            CHECK='Y'
            IF ((TYPE.EQ.'L').OR.(TYPE.EQ.'I').OR.(TYPE.EQ.'F').OR.     &
     &          (TYPE.EQ.'C')) CALL CHKTYP(WORK,DUMMY,TYPE,CHECK,FORM)
! invalid item here ---> stop counting items and say good bye ... !
            IF (CHECK.EQ.'N') GOTO 500
! 'unreadable' item ---> stop counting items and say good bye ... !
            IF ((CHECK.EQ.'U').AND.(TYPE.EQ.'I')) GOTO 500
            NITEMS=NITEMS+1
         ENDIF
         IF (IC.EQ.1) THEN
            CALL PARSE(WORK,BUFFER,NUMBER,'*',0)
            CHECK='Y'
            IF ((TYPE.EQ.'L').OR.(TYPE.EQ.'I').OR.(TYPE.EQ.'F').OR.     &
     &          (TYPE.EQ.'C')) CALL CHKTYP(NUMBER,DUMMY,TYPE,CHECK,FORM)
! invalid item here ---> stop counting items and say good bye ... !
            IF (CHECK.EQ.'N') GOTO 500
! 'unreadable' item ---> stop counting items and say good bye ... !
            IF ((CHECK.EQ.'U').AND.(TYPE.EQ.'I')) GOTO 500
            IF (EXPAND) CALL STRIP(NUMBER,L0,'A')
            CALL STRIP(BUFFER,L,'A')
! invalid item here ---> stop counting items and say good bye ... !
            IF (L.LT.1) GOTO 500
! 'multiplier'  MUST  be a positive integer number (strictly!) ...
            CALL CHKINT(BUFFER,DUMMY,CHECK,FORM)
! invalid item here ---> stop counting items and say good bye ... !
            IF (CHECK.NE.'Y') GOTO 500
! hopefully should work without error/end condition here ... ?



            DUMMY='('//FORM(1:253)//')'

            CALL STRIP(DUMMY,IMULT,'A')
            READ(BUFFER,DUMMY) IMULT
! invalid item here ---> stop counting items and say good bye ... !
            IF (IMULT.LE.0) GOTO 500
            IF (EXPAND) THEN
! invalid item here ---> stop counting items and say good bye ... !
               IF (L0.LT.1) GOTO 500
               IF (I.EQ.1) WORK=' '
               IF (I.GT.1) CALL SUBWRD(STRING,WORK,1,I-1)
               L=LENGTH(WORK)
               DO 200 J=1,IMULT
                  IF (I.EQ.1) WORK=NUMBER
                  IF (I.GT.1) WORK=WORK(1:L)//' '//NUMBER
                  L=L+L0+1
                  IF (I.EQ.1) L=L-1
                  IF ((L.GT.LW).OR.(L.GT.LS)) THEN
                     WRITE(*,'(A)') ' '
                     WRITE(*,'(A)') 'Error LEXLIB routine'//            &
     &                              ' ''NITEMS'': expansion fails,'//   &
     &                              ' insufficient CHARACTER length.'
! have to stop expansion ---> stop counting items and say good bye ... !
                     GOTO 500
                  ENDIF
  200          CONTINUE
               DO 300 J=I+1,N
                  CALL SUBWRD(STRING,BUFFER,J,1)
                  CALL STRIP(BUFFER,L0,'A')
                  IF (L0.LT.1) GOTO 500
                  WORK=WORK(1:L)//' '//BUFFER
                  L=L+L0+1
                  IF ((L.GT.LW).OR.(L.GT.LS)) THEN
                     WRITE(*,'(A)') ' '
                     WRITE(*,'(A)') 'Error LEXLIB routine'//            &
     &                              ' ''NITEMS'': expansion fails,'//   &
     &                              ' insufficient CHARACTER length.'
! have to stop expansion ---> stop counting items and say good bye ... !
                     GOTO 500
                  ENDIF
  300          CONTINUE
! sorry, do not know a simpler way to get new correct 'counters' here:
               STRING=WORK
               GOTO 100
            ELSE
               NITEMS=NITEMS+IMULT
            ENDIF
         ENDIF
  400 CONTINUE
  500 CONTINUE
      RETURN
      END


      SUBROUTINE SUBWRD(STRING,WORDS,IBEG,INUM)
      USE mod_comp
      IMPLICIT REAL(DP) (A-H,O-Z)
! Extracts specified words out of a string ...
      CHARACTER*(*) STRING,WORDS
      LOGICAL BLANK
      INTEGER NWORDS,L,LEN,LENGTH,I,IBEG,INUM,ISTART,ISTOP
      EXTERNAL LENGTH

      L=LENGTH(STRING)
      ISTART=0
      ISTOP=L
      NWORDS=0
      BLANK=.TRUE.
      DO 10 I=1,L
         IF (BLANK.AND.(STRING(I:I).NE.' ')) THEN
            BLANK=.FALSE.
            NWORDS=NWORDS+1
            IF (NWORDS.EQ.IBEG) ISTART=I
         ELSE IF ((.NOT.BLANK).AND.(STRING(I:I).EQ.' ')) THEN
            BLANK=.TRUE.
            IF (NWORDS.EQ.(IBEG+INUM-1)) ISTOP=I-1
         ENDIF
   10 CONTINUE
      IF ((ISTART.GT.0).AND.(INUM.GT.0)) THEN
         WORDS=STRING(ISTART:ISTOP)
         IF ((ISTOP-ISTART+1).GT.LEN(WORDS)) THEN
            WRITE(*,'(A)') ' '
            WRITE(*,'(A)') 'Warning LEXLIB routine ''SUBWRD'': '//      &
     &             'Output string will be truncated! The output is'
            WRITE(*,'(A,I5,A,I5,A)') 'a string of length ',             &
     &             ISTOP-ISTART+1,                                      &
     &             ' characters but ''WORDS'' can only hold ',          &
     &             LEN(WORDS),' characters.'
            WRITE(*,'(A)') 'Continuing execution ...'
            WRITE(*,'(A)') ' '
         ENDIF
      ELSE
         WORDS=' '
      ENDIF
      RETURN
      END


      SUBROUTINE STRIP(STRING,L,MODE)
      USE mod_comp
      IMPLICIT REAL(DP) (A-H,O-Z)
! Strips off blanks in STRING according to setting of MODE and returns
! the position L of the last non-blank character after all operations.
! MODE may be set to:
!   - 'L' remove all leading blanks only
!   - 'I' remove all blanks inside STRING, let leading blanks untouched
!   - 'S' merge all multiple blanks within STRING into one single blank
!         but leave all leading blanks untouched!
!   - 'B' remove leading blanks and merge all multiple blanks into one
!   - 'A' remove all (but really  all!) blanks
! all other settings lead to output L=0 (returns a blank string)!
      CHARACTER*(*) STRING
      CHARACTER*1   MODE
      INTEGER       L,LENGTH,L0,FIRST,POS,I
      EXTERNAL      LENGTH

      L0=LENGTH(STRING)

! Here stripping off all leading blanks:
      IF ((MODE.EQ.'L').OR.(MODE.EQ.'A').OR.(MODE.EQ.'B')) THEN
         FIRST=L0+1
         DO 10 I=1,L0
            IF (STRING(I:I).NE.' ') THEN
               FIRST=I
               GOTO 20
            ENDIF
  10     CONTINUE
  20     CONTINUE
         IF (FIRST.LE.L0) THEN
            STRING=STRING(FIRST:L0)
            L=L0
            IF (MODE.EQ.'L') L=LENGTH(STRING)
         ELSE
            STRING=' '
            L=0
         ENDIF
      END IF
! Here stripping off all blanks inside STRING except for leading blanks
      IF ((MODE.EQ.'I').OR.(MODE.EQ.'A')) THEN
         FIRST=L0+1
         DO 30 I=1,L0
            IF (STRING(I:I).NE.' ') THEN
               FIRST=I
               GOTO 40
            ENDIF
  30     CONTINUE
  40     CONTINUE
         IF (FIRST.LE.L0) THEN
            POS=FIRST+1
            DO 50 I=FIRST+1,L0
               IF (STRING(POS:POS).EQ.' ') THEN
                  IF (POS.LT.L0) THEN
                     STRING=STRING(1:POS-1)//STRING(POS+1:L0)
                  ELSE
                     STRING=STRING(1:POS-1)
                  ENDIF
               ELSE
                  POS=POS+1
               ENDIF
   50       CONTINUE
            L=LENGTH(STRING)
         ELSE
            STRING=' '
            L=0
         ENDIF
      ENDIF
! Here merging multiple blanks into a single blank (except leading blanks)
      IF ((MODE.EQ.'S').OR.(MODE.EQ.'B')) THEN
         FIRST=L0+1
         DO 60 I=1,L0
            IF (STRING(I:I).NE.' ') THEN
               FIRST=I
               GOTO 70
            ENDIF
  60     CONTINUE
  70     CONTINUE
         IF (FIRST.LE.L0) THEN
            POS=FIRST+1
            DO 80 I=FIRST+1,L0-1
               IF (STRING(POS:POS+1).EQ.'  ') THEN
                  IF (POS.LT.L0) THEN
                     STRING=STRING(1:POS-1)//STRING(POS+1:L0)
                  ELSE
                     STRING=STRING(1:POS-1)
                  ENDIF
               ELSE
                  POS=POS+1
               ENDIF
   80       CONTINUE
            L=LENGTH(STRING)
         ELSE
            STRING=' '
            L=0
         ENDIF
      ENDIF
! Invalid mode --> return blank string
      IF ((MODE.NE.'L').AND.(MODE.NE.'I').AND.(MODE.NE.'S')             &
     &                 .AND.(MODE.NE.'A').AND.(MODE.NE.'B')) THEN
         STRING=' '
         L=0
      ENDIF

      RETURN
      END


      SUBROUTINE UPPER(STRING)
      USE mod_comp
      IMPLICIT REAL(DP) (A-H,O-Z)
! uppercase all letters in STRING ...
      CHARACTER*(*) STRING
      CHARACTER*1   ALPHAL(26),ALPHAU(26)
      INTEGER       L,I,J,LENGTH
      EXTERNAL      LENGTH
      SAVE          ALPHAL,ALPHAU
      DATA ALPHAL /'a','b','c','d','e','f','g','h','i','j','k','l','m', &
     &             'n','o','p','q','r','s','t','u','v','w','x','y','z'/ 
      DATA ALPHAU /'A','B','C','D','E','F','G','H','I','J','K','L','M', &
     &             'N','O','P','Q','R','S','T','U','V','W','X','Y','Z'/

      L=LENGTH(STRING)
      DO 30 I=1,L
         DO 10 J=1,26
            IF (STRING(I:I).EQ.ALPHAL(J)) THEN
               STRING(I:I)=ALPHAU(J)
               GOTO 20
            ENDIF
   10    CONTINUE
   20    CONTINUE
   30 CONTINUE

      RETURN
      END


      SUBROUTINE LOWER(STRING)
      USE mod_comp
      IMPLICIT REAL(DP) (A-H,O-Z)
! lowercase all letters in STRING ...
      CHARACTER*(*) STRING
      CHARACTER*1   ALPHAL(26),ALPHAU(26)
      INTEGER       L,I,J,LENGTH
      EXTERNAL      LENGTH
      SAVE          ALPHAL,ALPHAU
      DATA ALPHAL /'a','b','c','d','e','f','g','h','i','j','k','l','m', &
     &             'n','o','p','q','r','s','t','u','v','w','x','y','z'/
      DATA ALPHAU /'A','B','C','D','E','F','G','H','I','J','K','L','M', &
     &             'N','O','P','Q','R','S','T','U','V','W','X','Y','Z'/

      L=LENGTH(STRING)
      DO 30 I=1,L
         DO 10 J=1,26
            IF (STRING(I:I).EQ.ALPHAU(J)) THEN
               STRING(I:I)=ALPHAL(J)
               GOTO 20
            ENDIF
   10    CONTINUE
   20    CONTINUE
   30 CONTINUE

      RETURN
      END


      SUBROUTINE PARSE(STRING,BEFORE,AFTER,PATTERN,MODE)
      USE mod_comp
      IMPLICIT REAL(DP) (A-H,O-Z)
! Parses STRING like REXX-interpreter with rule BEFORE 'PATTERN' AFTER;
! MODE determines which "length definition" shall be used for PATTERN
! (>=0: take length from LEN, means treat also trailing blanks or <0:
! take length from LENGTH, means ignore all trailing blanks ...).
      CHARACTER*(*) STRING,BEFORE,AFTER,PATTERN
      INTEGER       MODE,LENGTH,LS,LP,I
      EXTERNAL      LENGTH

      BEFORE=STRING
      AFTER=' '
      IF (MODE.GE.0) THEN
         LP=LEN(PATTERN)
      ELSE
! additional remark: blank string shall act like "single blank" ...
         LP=MAX(LENGTH(PATTERN),1)
      ENDIF
      LS=LENGTH(STRING)
      I=INDEX(STRING(1:LS),PATTERN(1:LP))
      IF (I.EQ.0) RETURN
      IF ((I.GT.1).AND.(I+LP.LE.LS)) THEN
         BEFORE=STRING(1:I-1)
         AFTER=STRING(I+LP:LS)
      ELSE IF ((I.EQ.1).AND.(I+LP.LE.LS)) THEN
         BEFORE=' '
         AFTER=STRING(I+LP:LS)
      ELSE IF ((I.EQ.1).AND.(I+LP.GT.LS)) THEN
         BEFORE=' '
         AFTER=' '
      ELSE IF (I+LP.GT.LS) THEN
         BEFORE=STRING(1:I-1)
         AFTER=' '
      ENDIF

      RETURN
      END


      FUNCTION NOCCUR(STRING,PATTERN,MODE)
      USE mod_comp
      IMPLICIT REAL(DP) (A-H,O-Z)
! Tells how often PATTERN occurs within STRING!
! MODE determines which "length definition" shall be used for PATTERN
! (>=0: take length from LEN, means treat also trailing blanks or <0:
! take length from LENGTH, means ignore all trailing blanks ...).
      CHARACTER*(*) STRING,PATTERN
      INTEGER       NOCCUR,LENGTH,LS,LP,I,LAST,MODE
      EXTERNAL      LENGTH

      NOCCUR=0
      LAST=1
      LS=LENGTH(STRING)
      IF (MODE.GE.0) THEN
         LP=LEN(PATTERN)
      ELSE
! again blank strings should be interpreted as 'single blank':
         LP=MAX(LENGTH(PATTERN),1)
      ENDIF
   10 CONTINUE
      I=INDEX(STRING(LAST:LS),PATTERN(1:LP))
      IF (I.EQ.0) RETURN
      NOCCUR=NOCCUR+1
      LAST=LAST+I+LP-1
      IF ((LAST+LP).GT.LS) RETURN
      GOTO 10

      RETURN
      END


      FUNCTION INDEXN(STRING,PATTERN,NTH)
      USE mod_comp
      IMPLICIT REAL(DP) (A-H,O-Z)
! Get the starting position of PATTERN within STRING for the NTHth
! occurence (it is some "generalized INDEX-function" ...):
      CHARACTER*(*) STRING,PATTERN
      INTEGER       INDEXN,NTH,LENGTH,LS,LP,I,LAST,NOCCUR
      EXTERNAL      LENGTH

      INDEXN=0
      NOCCUR=0
      LAST=1
      LS=LENGTH(STRING)
      IF (NTH.GE.0) THEN
         LP=LEN(PATTERN)
      ELSE
         LP=LENGTH(PATTERN)
      ENDIF
   10 CONTINUE
      I=INDEX(STRING(LAST:LS),PATTERN(1:LP))
      IF (I.EQ.0) RETURN
      NOCCUR=NOCCUR+1
      IF (NOCCUR.EQ.ABS(NTH)) INDEXN=LAST+I-1
      LAST=LAST+I+LP-1
      IF ((NOCCUR.EQ.ABS(NTH)).OR.((LAST+LP).GT.LS)) RETURN
      GOTO 10

      RETURN
      END


      SUBROUTINE REPLAC(STRING,OLDPAT,NEWPAT,L,MODE)
      USE mod_comp
      IMPLICIT REAL(DP) (A-H,O-Z)
! Replace some pattern OLDPAT by pattern NEWPAT within STRING according
! to the setting of MODE which may take following values:
!   -  0  global replacement (for all --- but really all --- occurences)
!   -  >0 replacement for the MODEth occurence only
!   -  <0 replacement for the first ABS(MODE) occurences
! on output L returns the position of the last non-blank character
! in STRING after all replacements, on input it controls which 'length'
! for OLDPAT/NEWPAT should be taken (L>0: length defined by LEN, means
! treat also trailing blanks or L<0: length defined by LENGTH, means all
! trailing blanks will be ignored ... . L=0 has the special meaning: ``do
! the same as for L>0 for OLDPAT, but the same as L<0 for NEWPAT'' (will
! become important for the special case when one wants exact replacement
! of OLDPAT including trailing blanks by a 'null string' ---> NEWPAT=' '
! and L=0 should be taken if one wishes to do that ...).
      CHARACTER*(*) STRING,OLDPAT,NEWPAT
      INTEGER      L,MODE,LENGTH,LOLD,LNEW,N,N1,N2,NOCCUR,I,INDEXN,L0,J
      INTEGER      IU,NXTFRU
      EXTERNAL     LENGTH,NOCCUR,INDEXN,NXTFRU

      L0=LEN(STRING)
      IF (L.GT.0) THEN
         J=1
         LOLD=LEN(OLDPAT)
         LNEW=LEN(NEWPAT)
      ELSE IF (L.EQ.0) THEN
         J=1
         LOLD=LEN(OLDPAT)
         LNEW=LENGTH(NEWPAT)
      ELSE
         J=-1
         LOLD=LENGTH(OLDPAT)
         LNEW=LENGTH(NEWPAT)
      ENDIF
! Do not replace 'null strings' (how to do??)!
      IF (LOLD.EQ.0) RETURN
      L=LENGTH(STRING)
! Nothing what could be replaced ...
      IF (L.EQ.0) RETURN
      N=NOCCUR(STRING,OLDPAT(1:LOLD),J)
      IF (N.EQ.0) RETURN
      IF (MODE.LT.0) THEN
         N1=1
         N2=ABS(MODE)
         IF (N2.GT.N) N2=N
      ELSE IF (MODE.GT.0) THEN
         IF (MODE.GT.N) RETURN
         N1=MODE
         N2=MODE
      ELSE
         N1=1
         N2=N
      ENDIF
      DO 10 N=N1,N2
         I=INDEXN(STRING,OLDPAT(1:LOLD),J*N1)
         IF ((LNEW.GT.0).AND.(LNEW.LE.LOLD)) THEN
            IF ((I.GT.1).AND.((I+LOLD).LE.L0)) THEN
               STRING=STRING(1:I-1)//NEWPAT(1:LNEW)//STRING(I+LOLD:L0)
            ELSE IF (I.EQ.1) THEN
               STRING=NEWPAT(1:LNEW)//STRING(1+LOLD:L0)
            ELSE IF ((I+LOLD).GT.L0) THEN
               STRING=STRING(1:I-1)//NEWPAT(1:LNEW)
            ENDIF
            L=L+LNEW-LOLD
         ELSE IF ((LNEW.GT.0).AND.(LNEW.GT.LOLD)) THEN
! Here we run into trouble due to the order in which the expression and
! the assignment are done: some intermediate partially replaced string
! will be used for the further evaluation of the expression!! So we have
! to take this into account very carefully (if LNEW>LOLD) ... !! Shit!!!
! We do not want to waste space for another temporary character variable
! so we do it the ugly way: via external I/O (on a scratch file). Sorry!
            IU=NXTFRU()
            OPEN(IU,STATUS='SCRATCH')
            IF ((I.GT.1).AND.((I+LOLD).LE.L0)) THEN
               WRITE(IU,'(A,A,A)')                                      &
     &                 STRING(1:I-1),NEWPAT(1:LNEW),STRING(I+LOLD:L0)
            ELSE IF (I.EQ.1) THEN
               WRITE(IU,'(A,A)') NEWPAT(1:LNEW),STRING(1+LOLD:L0)
            ELSE IF ((I+LOLD).GT.L0) THEN
               WRITE(IU,'(A,A)') STRING(1:I-1),NEWPAT(1:LNEW)
            ENDIF
            REWIND IU
            READ(IU,'(A)') STRING
            CLOSE(IU)
            L=L+LNEW-LOLD
            IF (L.GT.L0) THEN
               WRITE(*,'(A)') ' '
               WRITE(*,'(A)') 'Warning LEXLIB routine ''REPLAC'': '//   &
     &                'Output string will be truncated! New string is'
               WRITE(*,'(A)') 'longer than old string and the length'// &
     &                    ' of variable ''STRING'' is too short! Length'
               WRITE(*,'(A,I6,A,I5,A)') 'of new string is ',L,          &
     &                ' characters but ''STRING'' can only hold ',L0,   &
     &                ' characters.'
               WRITE(*,'(A)') 'Continuing execution ...'
               WRITE(*,'(A)') ' '
            ENDIF
         ELSE
! special code for 'null string' replacement ...
            IF ((I.GT.1).AND.((I+LOLD).LE.L0)) THEN
               STRING=STRING(1:I-1)//STRING(I+LOLD:L0)
            ELSE IF (I.EQ.1) THEN
               STRING=STRING(1+LOLD:L0)
            ELSE IF ((I+LOLD).GT.L0) THEN
               STRING=STRING(1:I-1)
            ENDIF
            L=L-LOLD
         ENDIF
   10 CONTINUE
      L=MIN(L,L0)

      RETURN
      END


      SUBROUTINE CHKTYP(STRING,WORK,MATCH,TYPE,FORM)
      USE mod_comp
      IMPLICIT REAL(DP) (A-H,O-Z)
! Try to check the validity of the data type of what is contained in
! STRING ... . The type to be tested must be given in MATCH and must be
!   - 'I'  string should contain a valid FORTRAN-Integer
!   - 'F'  string should contain a valid FORTRAN-Float (any format)
!   - 'C'  string should contain a FORTRAN-Complex number (any format)
!   - 'L'  string should contain a valid FORTRAN-Logical (any format)
!   - 'A'  string contains only alphanumeric characters [0-9,A-Z,a-z]
!   - 'U'  uppercase (alphanumeric) string
!   - 'l'  lowercase (alphanumeric) string
!   - 'H'  string would be a valid hexadecimal number (chars 0-9,A-F)
!   - 'O'  string would be a valid octal number (characters 0-7 only)
!   - 'B'  string would be a valid binary number (only 0's and 1's)
!   - 'N'  string is a 'null string' (empty string)
!   - any other (invalid) inputs will behave as 'test was negative'!
! The result is returned in TYPE ('N' means 'no' = false and all other
! means 'yes' = true [usually 'Y' returned, sometimes other values ...]
! In FORM a format string is returned being needed to read from STRING
! GENERAL WARNING: *all* blanks are ignored --> we test only 'one word'
      CHARACTER*(*) STRING,WORK,FORM
      CHARACTER*1   MATCH,TYPE,CH
      CHARACTER*15  FORM1,FORM2
      CHARACTER*255 PART1,PART2
      LOGICAL       LTEST,LPURE
      INTEGER       LENGTH,LEN,NOCCUR,L,I,J
      INTRINSIC     LEN
      EXTERNAL      LENGTH,NOCCUR

      TYPE='N'
      FORM=' '
      LTEST=.TRUE.
      LPURE=.TRUE.
      WORK=STRING
      CALL STRIP(WORK,L,'A')
      IF (MATCH.EQ.'I') THEN
         CALL CHKINT(STRING,WORK,TYPE,FORM)
         RETURN
      ELSE IF (MATCH.EQ.'F') THEN
         CALL CHKFLT(STRING,WORK,TYPE,FORM)
         RETURN
      ELSE IF (MATCH.EQ.'C') THEN
! a complex number must be of type '(' float/int ',' float/int ')'
         IF ((NOCCUR(WORK,',',0).NE.1).OR.(WORK(1:1).NE.'(')            &
     &          .OR.(WORK(L:L).NE.')').OR.(L.LE.4)) LTEST=.FALSE.
! all seems to be okay until here ...
         IF (LTEST) THEN
            CALL PARSE(WORK(2:L-1),PART1,PART2,',',0)
! part1 is float (or integer -- does not matter ...)?
            CALL CHKFLT(PART1,WORK,CH,FORM1)
            IF (CH.EQ.'N') LTEST=.FALSE.
            IF (LTEST) I=LENGTH(FORM1)
         ENDIF
! if part1 was okay, then we have still to test part2 ...
         IF (LTEST) THEN
! part2 is float (or integer -- does not matter ...)?
            CALL CHKFLT(PART2,WORK,CH,FORM2)
            IF (CH.EQ.'N') LTEST=.FALSE.
            IF (LTEST) J=LENGTH(FORM2)
         ENDIF
         IF (LTEST) THEN
            TYPE='Y'
! format for this input ...
            IF ((I.NE.0).AND.(J.NE.0).AND.(LEN(FORM).GE.(I+J+10)))      &
     &                 FORM='1X,'//FORM1(1:I)//',1X,'//FORM2(1:J)//',1X'
         ENDIF
         RETURN
      ELSE IF (MATCH.EQ.'L') THEN
! that is very easy (according to FORTRAN rules it must start with
! 'T', 'F' or '.T', '.F' --- and then may follow what you want:
         IF ((WORK(1:1).EQ.'T').OR.(WORK(1:1).EQ.'F').OR.               &
     &       (WORK(1:2).EQ.'.T').OR.(WORK(1:2).EQ.'.F')) TYPE='Y'
         IF ((TYPE.EQ.'Y').AND.(LEN(FORM).GE.6)) THEN
            WRITE(FORM,'(A1,I5)') 'L',L
            CALL STRIP(FORM,I,'A')
         ENDIF
         RETURN
      ELSE IF (MATCH.EQ.'A') THEN
! Check for 'empty string' ...
         IF (L.LE.0) THEN
            LTEST=.FALSE.
            GOTO 101
         ENDIF
! Case is not relevant ...
         CALL UPPER(WORK)
! Well, check ... :
         DO 100 I=1,L
            CH=WORK(I:I)
            IF ((CH.NE.'0').AND.(CH.NE.'1').AND.(CH.NE.'2').AND.        &
     &          (CH.NE.'3').AND.(CH.NE.'4').AND.(CH.NE.'5').AND.        &
     &          (CH.NE.'6').AND.(CH.NE.'7').AND.(CH.NE.'8').AND.        &
     &          (CH.NE.'9').AND.(CH.NE.'A').AND.(CH.NE.'B').AND.        &
     &          (CH.NE.'C').AND.(CH.NE.'D').AND.(CH.NE.'E').AND.        &
     &          (CH.NE.'F').AND.(CH.NE.'G').AND.(CH.NE.'H').AND.        &
     &          (CH.NE.'I').AND.(CH.NE.'J').AND.(CH.NE.'K').AND.        &
     &          (CH.NE.'L').AND.(CH.NE.'M').AND.(CH.NE.'N').AND.        &
     &          (CH.NE.'O').AND.(CH.NE.'P').AND.(CH.NE.'Q').AND.        &
     &          (CH.NE.'R').AND.(CH.NE.'S').AND.(CH.NE.'T').AND.        &
     &          (CH.NE.'U').AND.(CH.NE.'V').AND.(CH.NE.'W').AND.        &
     &          (CH.NE.'X').AND.(CH.NE.'Y').AND.(CH.NE.'Z'))            &
     &                                                    LTEST=.FALSE.
            IF (.NOT.LTEST) GOTO 101
            LPURE=LPURE.AND.(CH.NE.'0').AND.(CH.NE.'1').AND.            &
     &                      (CH.NE.'2').AND.(CH.NE.'3').AND.            &
     &                      (CH.NE.'4').AND.(CH.NE.'5').AND.            &
     &                      (CH.NE.'6').AND.(CH.NE.'7').AND.            &
     &                      (CH.NE.'8').AND.(CH.NE.'9')
  100    CONTINUE
  101    IF (LTEST) TYPE='Y'
! 'pure string' [no chars 0-9] shall return TYPE='P' instead of TYPE='Y'
         IF (LTEST.AND.LPURE) TYPE='P'
      ELSE IF (MATCH.EQ.'U') THEN
! Check for 'empty string' ...
         IF (L.LE.0) THEN
            LTEST=.FALSE.
            GOTO 201
         ENDIF
! Well, check ... :
         DO 200 I=1,L
            CH=WORK(I:I)
            IF ((CH.NE.'0').AND.(CH.NE.'1').AND.(CH.NE.'2').AND.        &
     &          (CH.NE.'3').AND.(CH.NE.'4').AND.(CH.NE.'5').AND.        &
     &          (CH.NE.'6').AND.(CH.NE.'7').AND.(CH.NE.'8').AND.        &
     &          (CH.NE.'9').AND.(CH.NE.'A').AND.(CH.NE.'B').AND.        &
     &          (CH.NE.'C').AND.(CH.NE.'D').AND.(CH.NE.'E').AND.        &
     &          (CH.NE.'F').AND.(CH.NE.'G').AND.(CH.NE.'H').AND.        &
     &          (CH.NE.'I').AND.(CH.NE.'J').AND.(CH.NE.'K').AND.        &
     &          (CH.NE.'L').AND.(CH.NE.'M').AND.(CH.NE.'N').AND.        &
     &          (CH.NE.'O').AND.(CH.NE.'P').AND.(CH.NE.'Q').AND.        &
     &          (CH.NE.'R').AND.(CH.NE.'S').AND.(CH.NE.'T').AND.        &
     &          (CH.NE.'U').AND.(CH.NE.'V').AND.(CH.NE.'W').AND.        &
     &          (CH.NE.'X').AND.(CH.NE.'Y').AND.(CH.NE.'Z'))            &
     &                                                    LTEST=.FALSE.
            IF (.NOT.LTEST) GOTO 201
            LPURE=LPURE.AND.(CH.NE.'0').AND.(CH.NE.'1').AND.            &
     &                      (CH.NE.'2').AND.(CH.NE.'3').AND.            &
     &                      (CH.NE.'4').AND.(CH.NE.'5').AND.            &
     &                      (CH.NE.'6').AND.(CH.NE.'7').AND.            &
     &                      (CH.NE.'8').AND.(CH.NE.'9')
  200    CONTINUE
  201    IF (LTEST) TYPE='Y'
! 'pure string' [no chars 0-9] shall return TYPE='P' instead of TYPE='Y'
         IF (LTEST.AND.LPURE) TYPE='P'
      ELSE IF (MATCH.EQ.'l') THEN
! Check for 'empty string' ...
         IF (L.LE.0) THEN
            LTEST=.FALSE.
            GOTO 301
         ENDIF
! Well, check ... :
         DO 300 I=1,L
            CH=WORK(I:I)
            IF ((CH.NE.'0').AND.(CH.NE.'1').AND.(CH.NE.'2').AND.        &
     &          (CH.NE.'3').AND.(CH.NE.'4').AND.(CH.NE.'5').AND.        &
     &          (CH.NE.'6').AND.(CH.NE.'7').AND.(CH.NE.'8').AND.        &
     &          (CH.NE.'9').AND.(CH.NE.'a').AND.(CH.NE.'b').AND.        &
     &          (CH.NE.'c').AND.(CH.NE.'d').AND.(CH.NE.'e').AND.        &
     &          (CH.NE.'f').AND.(CH.NE.'g').AND.(CH.NE.'h').AND.        &
     &          (CH.NE.'i').AND.(CH.NE.'j').AND.(CH.NE.'k').AND.        &
     &          (CH.NE.'l').AND.(CH.NE.'m').AND.(CH.NE.'n').AND.        &
     &          (CH.NE.'o').AND.(CH.NE.'p').AND.(CH.NE.'q').AND.        &
     &          (CH.NE.'r').AND.(CH.NE.'s').AND.(CH.NE.'t').AND.        &
     &          (CH.NE.'u').AND.(CH.NE.'v').AND.(CH.NE.'w').AND.        &
     &          (CH.NE.'x').AND.(CH.NE.'y').AND.(CH.NE.'z'))            &
     &                                                    LTEST=.FALSE.
            IF (.NOT.LTEST) GOTO 301
            LPURE=LPURE.AND.(CH.NE.'0').AND.(CH.NE.'1').AND.            &
     &                      (CH.NE.'2').AND.(CH.NE.'3').AND.            &
     &                      (CH.NE.'4').AND.(CH.NE.'5').AND.            &
     &                      (CH.NE.'6').AND.(CH.NE.'7').AND.            &
     &                      (CH.NE.'8').AND.(CH.NE.'9')
  300    CONTINUE
  301    IF (LTEST) TYPE='Y'
! 'pure string' [no chars 0-9] shall return TYPE='P' instead of TYPE='Y'
         IF (LTEST.AND.LPURE) TYPE='P'
      ELSE IF (MATCH.EQ.'H') THEN
! Check for 'empty string' ...
         IF (L.LE.0) THEN
            LTEST=.FALSE.
            GOTO 401
         ENDIF
! Case shall not be of interest ...
         CALL UPPER(WORK)
! Well, check ... :
         DO 400 I=1,L
            CH=WORK(I:I)
            IF ((CH.NE.'0').AND.(CH.NE.'1').AND.(CH.NE.'2').AND.        &
     &          (CH.NE.'3').AND.(CH.NE.'4').AND.(CH.NE.'5').AND.        &
     &          (CH.NE.'6').AND.(CH.NE.'7').AND.(CH.NE.'8').AND.        &
     &          (CH.NE.'9').AND.(CH.NE.'A').AND.(CH.NE.'B').AND.        &
     &          (CH.NE.'C').AND.(CH.NE.'D').AND.(CH.NE.'E').AND.        &
     &          (CH.NE.'F')) LTEST=.FALSE.
            IF (.NOT.LTEST) GOTO 401
  400    CONTINUE
  401    IF (LTEST) TYPE='Y'
      ELSE IF (MATCH.EQ.'O') THEN
! Check for 'empty string' ...
         IF (L.LE.0) THEN
            LTEST=.FALSE.
            GOTO 501
         ENDIF
! Well, check ... :
         DO 500 I=1,L
            CH=WORK(I:I)
            IF ((CH.NE.'0').AND.(CH.NE.'1').AND.(CH.NE.'2').AND.        &
     &          (CH.NE.'3').AND.(CH.NE.'4').AND.(CH.NE.'5').AND.        &
     &          (CH.NE.'6').AND.(CH.NE.'7')) LTEST=.FALSE.
            IF (.NOT.LTEST) GOTO 501
  500    CONTINUE
  501    IF (LTEST) TYPE='Y'
      ELSE IF (MATCH.EQ.'B') THEN
! Check for 'empty string' ...
         IF (L.LE.0) THEN
            LTEST=.FALSE.
            GOTO 601
         ENDIF
! Well, check ... :
         DO 600 I=1,L
            CH=WORK(I:I)
            IF ((CH.NE.'0').AND.(CH.NE.'1')) LTEST=.FALSE.
            IF (.NOT.LTEST) GOTO 601
  600    CONTINUE
  601    IF (LTEST) TYPE='Y'
      ELSE IF (MATCH.EQ.'N') THEN
! that is too easy ...
         IF (L.LE.0) TYPE='Y'
      ENDIF
! all what ends up here were checks for things which could only be read
! as strings (A-format for FORTRAN read ...), give here the format ...
      IF ((TYPE.NE.'N').AND.(LEN(FORM).GE.6).AND.(MATCH.NE.'N')) THEN
         WRITE(FORM,'(A1,I5)') 'A',L
         CALL STRIP(FORM,I,'A')
      ENDIF

      RETURN
      END


      SUBROUTINE DATTYP(STRING,WORK,TYPE,FORM)
      USE mod_comp
      IMPLICIT REAL(DP) (A-H,O-Z)
! Try to find out the data type of what is contained in STRING ... .
! Of course we do not want to distinguish too may special cases here!
! The result will be returned in TYPE which could take the values:
!   - 'G'  'general string' containing any arbitrary things ...
!   - 'A'  alphanumeric string (contains only characters [0-9,A-Z,a-z])
!   - 'F'  valid floating point number
!   - 'C'  valid complex number
!   - 'I'  valid integer number (would also match type 'F')
!   - 'L'  valid logical value (would also match type 'A')
!   - 'N'  'null string' (empty string)
! Additionally we try to return an appropriate format string in FORM
! GENERAL WARNING: *all* blanks are ignored --> we test only 'one word'
! Let us start ...
      CHARACTER*(*) STRING,WORK,FORM
      CHARACTER*1   TYPE,MATCH
      INTEGER       LENGTH,L,LEN
      INTRINSIC     LEN
      EXTERNAL      LENGTH

! Now check possible type for possible type -- beginning with the most
! special possibilities and ending up with more and more general types:
      TYPE='N'
! 'null string' ...
      CALL CHKTYP(STRING,WORK,TYPE,MATCH,FORM)
      IF (MATCH.NE.'N') RETURN
      TYPE='L'
! logical value ... (string beginning with 'F','T' or '.F','.T')
      CALL CHKTYP(STRING,WORK,TYPE,MATCH,FORM)
      IF (MATCH.NE.'N') RETURN
      TYPE='I'
! integer ... (strict form only ...)
      CALL CHKTYP(STRING,WORK,TYPE,MATCH,FORM)
      IF (MATCH.EQ.'Y') RETURN
      TYPE='C'
! complex number ...
      CALL CHKTYP(STRING,WORK,TYPE,MATCH,FORM)
      IF (MATCH.NE.'N') RETURN
      TYPE='F'
! floating point number ...
      CALL CHKTYP(STRING,WORK,TYPE,MATCH,FORM)
      IF (MATCH.NE.'N') RETURN
      TYPE='A'
! alphanumeric string ...
      CALL CHKTYP(STRING,WORK,TYPE,MATCH,FORM)
      IF (MATCH.NE.'N') RETURN
      TYPE='G'
! okay, arriving here means 'general string' ...
      IF (LEN(FORM).GE.6) THEN
         WRITE(FORM,'(A1,I5)') 'A',LENGTH(STRING)
         CALL STRIP(FORM,L,'A')
      ENDIF

      RETURN
      END


      SUBROUTINE CHKINT(STRING,WORK,TYPE,FORM)
      USE mod_comp
      IMPLICIT REAL(DP) (A-H,O-Z)
! Supplementary routine for checking validity of type INTEGER, the
! answer is returned in TYPE ('Y' or 'N'), FORM is a format string
! for this integer number (if it is one ...) for FORTRAN-reading ...
      CHARACTER*(*) STRING,WORK,FORM
      CHARACTER*1   TYPE,CH
      CHARACTER*16  FTEST,FTEMP
      CHARACTER*255 PART1,PART2
      LOGICAL       LTEST
      INTEGER       L,I,J,K,LENGTH,NOCCUR,INDEXN,LEN
      REAL          RTEST,MAXINT,TINY
      INTRINSIC     LEN
      EXTERNAL      LENGTH,NOCCUR,INDEXN
! TINY is the machine tolerance and MAXINT the largest possible integer
! number (here for 32-bit, generally 2**(bit-1) - 1) --> maybe customize
      PARAMETER (TINY=1.E-10_DP,MAXINT=2147483647.E0_DP)

! default values ...
      TYPE='N'
      FORM=' '
      LTEST=.TRUE.
      WORK=STRING
      CALL STRIP(WORK,L,'A')
! First character may be a sign (+/-) or a number ...
      CH=WORK(1:1)
! if it is a sign something non-blank MUST follow ...
      IF (((CH.EQ.'+').OR.(CH.EQ.'-')).AND.(L.LE.1)) LTEST=.FALSE.
      IF ((CH.NE.'+').AND.(CH.NE.'-').AND.(CH.NE.'0').AND.              &
     &    (CH.NE.'1').AND.(CH.NE.'2').AND.(CH.NE.'3').AND.              &
     &    (CH.NE.'4').AND.(CH.NE.'5').AND.(CH.NE.'6').AND.              &
     &    (CH.NE.'7').AND.(CH.NE.'8').AND.(CH.NE.'9')) LTEST=.FALSE.
      IF (.NOT.LTEST) GOTO 101
      DO 100 I=2,L
! Now only numbers may follow ...
         CH=WORK(I:I)
         IF ((CH.NE.'0').AND.(CH.NE.'1').AND.(CH.NE.'2').AND.           &
     &       (CH.NE.'3').AND.(CH.NE.'4').AND.(CH.NE.'5').AND.           &
     &       (CH.NE.'6').AND.(CH.NE.'7').AND.(CH.NE.'8').AND.           &
     &       (CH.NE.'9')) LTEST=.FALSE.
         IF (.NOT.LTEST) GOTO 101
  100 CONTINUE
  101 CONTINUE
! now last check! is it within the allowed range for integer numbers???
      IF (LTEST) THEN
         WRITE(FTEST,'(A1,I5,A2)') 'F',L,'.0'
         FTEMP='('//FTEST(1:14)//')'
         CALL STRIP(FTEMP,I,'A')
         READ(WORK,FTEMP) RTEST
! sorry, number too large!
         IF (ABS(RTEST).GT.MAXINT) LTEST=.FALSE.
      ENDIF
! test was successful ...
      IF (LTEST) THEN
         TYPE='Y'
         IF (LEN(FORM).GE.6) THEN
            WRITE(FORM,'(A1,I5)') 'I',L
            CALL STRIP(FORM,I,'A')
         ENDIF
      ELSE
! maybe (and that is also some special case we might allow) it is a
! floating point number having an integer value, check it and if yes
! return TYPE='F' or 'E' (for 'F'/'E'-format ...) instead of TYPE='Y'
! but: allow only some maximum number MAXINT ("valid integer range")
         CALL CHKFLT(STRING,WORK,CH,FTEST)
         IF (CH.NE.'N') THEN
! all right it is a simple floating point number ...
            FTEMP='('//FTEST(1:14)//')'
            CALL STRIP(FTEMP,I,'A')
            READ(STRING,FTEMP) RTEST
            RTEST=ABS(RTEST)
            IF (CH.EQ.'Y') THEN
! all right it is a simple floating point number and not too large?
               I=LENGTH(FTEST)
               IF ((FTEST(1:1).EQ.'F').AND.(FTEST(I-1:I).EQ.'.0').AND.  &
     &                                         (RTEST.LT.MAXINT)) THEN
! it must be Fxxx.0-format, then it is okay and quite simple ...
                  TYPE='F'
                  IF (LEN(FORM).GE.9) THEN
                     WRITE(FORM,'(A1,I5,A3)') 'I',I-1,',1X'
                     CALL STRIP(FORM,I,'A')
                  ENDIF
               ENDIF
            ELSE
! E-format ...
               IF (((ABS(RTEST-NINT(RTEST))/RTEST).LT.TINY).AND.        &
     &                                         (RTEST.LT.MAXINT)) THEN
! bingo? -- it represents some (not too large) integer number ...
                  WORK=STRING
                  CALL STRIP(WORK,L,'A')
                  IF ((CH.EQ.'E').OR.(CH.EQ.'D').OR.(CH.EQ.'Q'))        &
     &                               CALL PARSE(WORK,PART1,PART2,CH,0)
                  IF (CH.EQ.'S') THEN
                     I=NOCCUR(WORK,'-',0)
                     J=NOCCUR(WORK,'+',0)
                     K=MAX(INDEXN(WORK,'-',I),INDEXN(WORK,'+',J))
                     PART1=WORK(1:K-1)
                     PART2=WORK(K:L)
                  ENDIF
                  CALL STRIP(PART1,I,'A')
                  CALL STRIP(PART2,J,'A')
! now some problem ... --> might get trouble finding some I-format; can
! only do it with constructs like yyyyyy(.)E(+/-)0, rest is impossible
! without changing string ... (so rest returns FORM=' ' and TYPE='U' !)
                  READ(PART2,'(I6)') K
                  IF (K.EQ.0) THEN
! exponent is zero, if it is now a Exxx.0-format it is okay and simple
                     K=LENGTH(FTEST)
                     IF (FTEST(K-1:K).EQ.'.0') THEN
                        TYPE='E'
                        IF (LEN(FORM).GE.13) THEN
                           I=I-NOCCUR(PART1,'.',0)
                           IF (CH.NE.'S') J=J+1
                           J=J+NOCCUR(PART1,'.',0)
                           WRITE(FORM,'(A1,I5,A1,I5,A1)')'I',I,',',J,'X'
                           CALL STRIP(FORM,J,'A')
                        ENDIF
                     ELSE
! sorry no I-format available for this number ... (should never occur?)
                        TYPE='U'
                     ENDIF
                  ELSE
! sorry no I-format available for this number ...
                     TYPE='U'
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDIF

      RETURN
      END


      SUBROUTINE CHKFLT(STRING,WORK,TYPE,FORM)
      USE mod_comp
      IMPLICIT REAL(DP) (A-H,O-Z)
! Supplementary routine for checking validity of floating point number,
! the answer is returned in TYPE ('Y' or 'N'), FORM is a format string
! for this floating point number (if it is one ...) for a FORTRAN-read
      CHARACTER*(*) STRING,WORK,FORM
      CHARACTER*1   TYPE,CH
      LOGICAL       LTEST,LDOT
      INTEGER       L,I,J,LEN
      INTRINSIC     LEN

! default values ...
      TYPE='N'
      FORM=' '
      LTEST=.TRUE.
      WORK=STRING
      CALL STRIP(WORK,L,'A')
      LDOT=.FALSE.
! First character may be a sign (+/-) a dot (.) or a number ...
      CH=WORK(1:1)
! if it is a sign or a dot something non-blank MUST follow ...
      IF (((CH.EQ.'+').OR.(CH.EQ.'-').OR.(CH.EQ.'.'))                   &
     &                              .AND.(L.LE.1)) LTEST=.FALSE.
      IF ((CH.NE.'+').AND.(CH.NE.'-').AND.(CH.NE.'0').AND.              &
     &    (CH.NE.'1').AND.(CH.NE.'2').AND.(CH.NE.'3').AND.              &
     &    (CH.NE.'4').AND.(CH.NE.'5').AND.(CH.NE.'6').AND.              &
     &    (CH.NE.'7').AND.(CH.NE.'8').AND.(CH.NE.'9').AND.              &
     &    (CH.NE.'.')) LTEST=.FALSE.
! only one single dot may be there, so remind if we had some already
      IF (CH.EQ.'.') THEN
         LDOT=.TRUE.
! the format would be the following (if the rest is correct ...)
         IF (LEN(FORM).GE.12) THEN
            WRITE(FORM,'(A1,I5,A1,I5)') 'F',L,'.',L-1
            CALL STRIP(FORM,J,'A')
         ENDIF
      ENDIF
      IF (.NOT.LTEST) GOTO 101
      DO 100 I=2,L
! Now only numbers may follow and somewhere a dot ...
         CH=WORK(I:I)
         IF ((CH.NE.'0').AND.(CH.NE.'1').AND.(CH.NE.'2').AND.           &
     &       (CH.NE.'3').AND.(CH.NE.'4').AND.(CH.NE.'5').AND.           &
     &       (CH.NE.'6').AND.(CH.NE.'7').AND.(CH.NE.'8').AND.           &
     &       (CH.NE.'9').AND.(CH.NE.'.')) LTEST=.FALSE.
! only one single dot may be there and remind if we had some already
         IF ((CH.EQ.'.').AND.LDOT) LTEST=.FALSE.
         IF (CH.EQ.'.') THEN
            LDOT=.TRUE.
! the format would be the following (if the rest is correct ...)
            IF (LEN(FORM).GE.12) THEN
               WRITE(FORM,'(A1,I5,A1,I5)') 'F',L,'.',L-I
               CALL STRIP(FORM,J,'A')
            ENDIF
         ENDIF
         IF (.NOT.LTEST) GOTO 101
  100 CONTINUE
! okay, LTEST should be true and LDOT too (otherwise it is an INTEGER!)
  101 IF (LTEST.AND.LDOT) TYPE='Y'
! if LTEST is true but LDOT is false we got an integer -- of course an
! integer is also valid for assignments to real values, so in this case
! we will not answer 'N' (no), but also not 'Y' (yes) --> answer is 'I'
      IF (LTEST.AND.(.NOT.LDOT)) THEN
         TYPE='I'
         IF (LEN(FORM).GE.12) THEN
            WRITE(FORM,'(A1,I5,A1,I5)') 'F',L,'.',0
            CALL STRIP(FORM,J,'A')
         ENDIF
      ENDIF
      IF (TYPE.EQ.'N') THEN
! if it is not a floating point number in direct F-format, it could be
! a number in 'E'-,'D'- or 'Q'-format ... --> check this here, it is of
! course also a correct and valid definition of a floating point number
         CALL CHKEDQ(STRING,WORK,'E',CH,FORM)
! return type='E' instead of type='Y'
         IF (CH.EQ.'Y') TYPE='E'
! or was it this ugly strange format without any 'E'/'D'/'Q'??
         IF (CH.EQ.'S') TYPE='S'
         IF (CH.EQ.'N') THEN
            CALL CHKEDQ(STRING,WORK,'D',CH,FORM)
! return type='D' instead of type='Y'
            IF (CH.EQ.'Y') TYPE='D'
            IF (CH.EQ.'N') THEN
               CALL CHKEDQ(STRING,WORK,'Q',CH,FORM)
! return type='Q' instead of type='Y'
               IF (CH.EQ.'Y') TYPE='Q'
            ENDIF
         ENDIF
      ENDIF

      RETURN
      END


      SUBROUTINE CHKEDQ(STRING,WORK,MATCH,TYPE,FORM)
      USE mod_comp
      IMPLICIT REAL(DP) (A-H,O-Z)
! Supplementary routine for checking validity of floating point numbers
! in the 'E', 'D' or 'Q'-format (if allowed/available ...) ...; the
! answer is returned in TYPE ('Y' or 'N'), FORM is a format string for
! this floating point number (if it is one ...) for FORTRAN-reading
      CHARACTER*(*) STRING,WORK,FORM
      CHARACTER*1   TYPE,CH,MATCH
      CHARACTER*255 PART1,PART2
      LOGICAL       LTEST,LDOT,ALLOWQ
      INTEGER       L,I,J,K,LEN,NOCCUR,LENGTH,INDEX,INDEXN
      INTRINSIC     LEN,INDEX,MAX
      EXTERNAL      NOCCUR,LENGTH,INDEXN
! Q-format (quadruple precision) is only supported on few machines ...
      PARAMETER (ALLOWQ=.TRUE.)

! default values ...
      TYPE='N'
      FORM=' '
! is 'Q' allowed/possile on this machine?
      IF ((MATCH.EQ.'Q').AND.(.NOT.ALLOWQ)) RETURN
! wrong matching type ...
      IF ((MATCH.NE.'E').AND.(MATCH.NE.'D').AND.(MATCH.NE.'Q')) RETURN
      LTEST=.TRUE.
      LDOT=.FALSE.
      K=0
      WORK=STRING
      CALL STRIP(WORK,L,'A')
! a float in E-format might contain one single 'E' or 'e', what is in
! front should be a float or an integer, what is after it an integer:
! (for D-format or Q-format just replace 'E' by 'D' or 'Q' ...). Of
! course there is also another valid (but more strange format without
! any 'E','D' or 'Q' -- just a signed exponent with some float in front;
! support also this valid format by searching such constructs ...
      CALL UPPER(WORK)
      IF (NOCCUR(WORK,MATCH,0).GT.1) LTEST=.FALSE.
! 'usual' format with 'E','D' or 'Q' here ...
      IF (LTEST.AND.(NOCCUR(WORK,MATCH,0).EQ.1)) THEN
! parse the result string ...
         CALL PARSE(WORK,PART1,PART2,MATCH,0)
         CALL STRIP(PART1,J,'A')
      ELSE IF (LTEST) THEN
! here the more strange (but also valid) format ..., it must at least
! contain one '+' or one '-' and at maximum two signs at all (the last
! is then relevant ...), so check all this here first ...
         I=NOCCUR(WORK,'-',0)
         J=NOCCUR(WORK,'+',0)
         IF (((I+J).NE.1).AND.((I+J).NE.2)) LTEST=.FALSE.
         IF (LTEST) THEN
! still all okay --> get last sign ...
            K=MAX(INDEXN(WORK,'-',I),INDEXN(WORK,'+',J))
! can not be possible ...
            IF ((K.LE.1).OR.(K.EQ.L)) LTEST=.FALSE.
            IF (LTEST) THEN
! well, seems still all to be okay -- parse the string ...
               PART1=WORK(1:K-1)
               CALL STRIP(PART1,J,'A')
               PART2=WORK(K:L)
            ENDIF
         ENDIF
      ENDIF
! if all is okay up to now do the rest ...
      IF (LTEST) THEN
! first part should be a floating point number ...
         CH=PART1(1:1)
! if it is a sign or a dot something non-blank MUST follow ...
         IF (((CH.EQ.'+').OR.(CH.EQ.'-').OR.(CH.EQ.'.'))                &
     &                                 .AND.(L.LE.1)) LTEST=.FALSE.
         IF ((CH.NE.'+').AND.(CH.NE.'-').AND.(CH.NE.'0').AND.           &
     &       (CH.NE.'1').AND.(CH.NE.'2').AND.(CH.NE.'3').AND.           &
     &       (CH.NE.'4').AND.(CH.NE.'5').AND.(CH.NE.'6').AND.           &
     &       (CH.NE.'7').AND.(CH.NE.'8').AND.(CH.NE.'9').AND.           &
     &       (CH.NE.'.')) LTEST=.FALSE.
! only one single dot may be there, so remind if we had some already
         IF (CH.EQ.'.') LDOT=.TRUE.
         IF (.NOT.LTEST) GOTO 101
         DO 100 I=2,J
! Now only numbers may follow and somewhere a dot ...
            CH=PART1(I:I)
            IF ((CH.NE.'0').AND.(CH.NE.'1').AND.(CH.NE.'2').AND.        &
     &          (CH.NE.'3').AND.(CH.NE.'4').AND.(CH.NE.'5').AND.        &
     &          (CH.NE.'6').AND.(CH.NE.'7').AND.(CH.NE.'8').AND.        &
     &          (CH.NE.'9').AND.(CH.NE.'.')) LTEST=.FALSE.
! only one single dot may be there and remind if we had some already
            IF ((CH.EQ.'.').AND.LDOT) LTEST=.FALSE.
            IF (CH.EQ.'.') LDOT=.TRUE.
            IF (.NOT.LTEST) GOTO 101
  100    CONTINUE
  101    CONTINUE
      ENDIF
      IF (LTEST) THEN
! if all is okay until here check finally the exponent (must be integer)
         CALL STRIP(PART2,J,'A')
! First character may be a sign (+/-) or a number ...
         CH=PART2(1:1)
! if it is a sign something non-blank MUST follow ...
         IF (((CH.EQ.'+').OR.(CH.EQ.'-')).AND.(L.LE.1)) LTEST=.FALSE.
         IF ((CH.NE.'+').AND.(CH.NE.'-').AND.(CH.NE.'0').AND.           &
     &       (CH.NE.'1').AND.(CH.NE.'2').AND.(CH.NE.'3').AND.           &
     &       (CH.NE.'4').AND.(CH.NE.'5').AND.(CH.NE.'6').AND.           &
     &       (CH.NE.'7').AND.(CH.NE.'8').AND.(CH.NE.'9')) LTEST=.FALSE.
         IF (.NOT.LTEST) GOTO 201
         DO 200 I=2,J
! Now only numbers may follow ...
            CH=PART2(I:I)
            IF ((CH.NE.'0').AND.(CH.NE.'1').AND.(CH.NE.'2').AND.        &
     &          (CH.NE.'3').AND.(CH.NE.'4').AND.(CH.NE.'5').AND.        &
     &          (CH.NE.'6').AND.(CH.NE.'7').AND.(CH.NE.'8').AND.        &
     &          (CH.NE.'9')) LTEST=.FALSE.
            IF (.NOT.LTEST) GOTO 201
  200    CONTINUE
  201    CONTINUE
      ENDIF
! if all was okay get the format (use always E-format ...):
      IF (LTEST) THEN
         TYPE='Y'
! 'strange format' shall return TYPE='S' instead of TYPE='Y' ...
         IF (K.NE.0) TYPE='S'
         IF (LEN(FORM).GE.12) THEN
! 'field length' in format is L, what is number of significant digits?
            I=0
            IF (LDOT) I=LENGTH(PART1)-INDEX(PART1,'.')
            WRITE(FORM,'(A1,I5,A1,I5)') 'E',L,'.',I
            CALL STRIP(FORM,I,'A')
         ENDIF
      ENDIF

      RETURN
      END

      function errfc(x)
      USE mod_comp
      IMPLICIT REAL(DP) (A-H,O-Z)
      z=dabs(x)
      t=1d0/(1d0+0.5d0*z)
      errfc=t*exp(-z*z-1.26551223d0+t*(1.00002368d0+t*(.37409196d0+     &
     &  t*(.09678418d0+t*(-.18628806d0+t*(.27886807d0+t*(-1.13520398d0+ &
     &   t*(1.48851587d0+t*(-.82215223d0+t*.17087277d0)))))))))
      if(x.lt.0d0) errfc=2d0-errfc

      return
      end

      function errf(x)
      USE mod_comp
      IMPLICIT REAL(DP) (A-H,O-Z)
      errf=1._DP-errfc(x)
      end function

      INTEGER FUNCTION NXTFRU()
      USE mod_comp
      IMPLICIT REAL(DP) (A-H,O-Z)
! Find the next free unit number ...
      LOGICAL OCCUP
      INTEGER I
      NXTFRU=-1
! Usually the standard FORTRAN range for unit numbers is 0...99;
! on some systems also unit numbers beyond 99 might be allowed:
! if this is the case and you want to make use of it change it ...
      DO 100 I=0,99
! Units 0,5,6 are usually reserved for stderr,stdin,stdout ...
! If your system uses other standard I/O-units change it ... !
         IF ((I.LE.0).OR.(I.EQ.5).OR.(I.EQ.6)) GOTO 100
         INQUIRE(UNIT=I,OPENED=OCCUP)
         IF (.NOT.OCCUP) THEN
            NXTFRU=I
            RETURN
         END IF
  100 CONTINUE
      RETURN
      END
