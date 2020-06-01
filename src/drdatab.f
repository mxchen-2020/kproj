      ! VASP subroutine
!***********************************************************************
!                                                                      *
      SUBROUTINE RDATAB(LOPEN,FNAM,IU,WHAT,PCHAR,CCHAR,SCHAR,TYPE, &
     &                     INTRES,FLTRES,CMPRES,LOGRES,STR,N,NMAX,IERR)
      use mod_comp
      implicit real(DP) (A-H,O-Z)
!                                                                      *
!  This routine extracts data from some file (name is stored in FNAM). *
!  The format is totally free, i.e. there may be as much blanks as you *
!  like and the data may be located everywhere. The exact rules are:   *
!                                                                      *
!     -- no line/definition may be longer than 500   characters        *
!     -- no 'keyword' may be longer than 255 characters and the        *
!        'keyword' is the first what appears in each definition!       *
!     -- PCHAR seperates the 'keyword' to be searched (WHAT) and the   *
!        data list (usual recommendation is to use '=' or ':' ...)     *
!     -- all after the character specified in CCHAR (recommendation    *
!        is to use either '#' or '!' ...) is treated as a comment      *
!     -- the character given in SCHAR shall seperate all definitions   *
!        given in one line (so one may place more than one definition  *
!        within one single line in order have the possibility to group *
!        things easy together which have some common meaning/function; *
!        recommendation is to use ';' or ',' ...)                      *
!     -- different data within one list must be seperated by blanks    *
!     -- no intermixing of different data types is allowed in a list!! *
!     -- in string data all blanks are taken 'as is' (from character   *
!        PCHAR up to the end of line or the next separator or comment) *
!     -- else multiple blanks act like one blank within a list of data *
!        or like 'no blank' if inserted between the 'keyword' and the  *
!        character PCHAR or between character PCHAR and the data list  *
!     -- because the characters given in PCHAR, CCHAR, SCHAR have now  *
!        already some special meaning one needs some extra way to make *
!        them accessible in string data: if one wants to use them in   *
!        any input strings one has to precede this character with \    *
!        (similar to UNIX-rules ...), because \ has also a special     *
!        meaning \ must be written as \\ within string data ...        *
!     -- data lists can be continued over several lines (with the only *
!        restriction that the total input may not be longer than 32767 *
!        characters ...). Continuations must be marked by some \ at    *
!        the end of the previous line (but really at the end of the    *
!        line!! -- so lines to be continued may not contain comments!) *
!     -- the naming of the 'keywords' shall be case-insensitive (this  *
!        would mean: 'TOKEN1', 'token1' and 'Token1' define in fact    *
!        all the same keyword !!) -- internally uppercase will be used *
!     -- for multiple entries (identical 'keywords') always the first  *
!        occurence defines what is extracted, all other entries will   *
!        be ignored!                                                   *
!                                                                      *
!  Maybe this routine is not yet perfect, but well usuable ...         *
!                                                                      *
!  There may be specified five types (given in variable TYPE) for the  *
!  data to be extracted:                                               *
!                                                                      *
!     -- 'S'    read some string       --> result goes to variable STR *
!     -- 'I'    read integer data list --> result goes to array INTRES *
!     -- 'F'    read real    data list --> result goes to array FLTRES *
!     -- 'C'    read complex data list --> result goes to array CMPRES *
!     -- 'L'    read logical data list --> result goes to array LOGRES *
!                                                                      *
!  for TYPE=I,F,C and L the variable N returns the number of elements  *
!  found in the data list, if an error occurred (e.g. wrong format) N  *
!  will be zero (and IERR non-zero). For TYPE=S variable N returns the *
!  position of the last non-blank character in variable STR ... .      *
!                                                                      *
!  IERR returns an error code if any problem occured (normal: IERR=0)  *
!     -- IERR=1   no free I/O-unit found to open file                  *
!     -- IERR=2   OPEN error (e.g. file not found, ...)                *
!     -- IERR=3   'keyword' (WHAT) not found on specified file         *
!     -- IERR=4   invalid data type (TYPE)                             *
!     -- IERR=5   error reading/parsing data list (check format!)      *
!     -- IERR=6   cannot open scratch file for conversion of data      *
!                                                                      *
!***********************************************************************

      CHARACTER*(*)   FNAM,WHAT,STR
      CHARACTER*1     PCHAR,CCHAR,SCHAR,TYPE,BS
      CHARACTER*2     HIDDEN
      CHARACTER*255   KEY
      CHARACTER*255   BUFLIN
      CHARACTER*255   WORK,WPARSE
      INTEGER         N,IERR,INTRES(*),I,IU,L,NMAX,IKEY,LKEY,INOSL
      REAL(DP)            FLTRES(*)
      COMPLEX(DPC)         CMPRES(*)
      LOGICAL         LOGRES(*),FOUND,CONT,YINTST,LOPEN
      INTEGER         LENGTH,NWORDS,NXTFRU,IUSCR,LMAX,NREAD,NITEMS
! Following parameter should probably be customized if porting this
! program from one machine to another machine: Usually in FORTRAN77
! it is  not   allowed to use the '*'-format for internal reads or
! writes, but we try to use it (if possible!). If your OS/compiler
! allows internal read/writes with the '*'-format (like AIX on RS6000)
! then set YINTST=.TRUE., otherwise you must use YINTST=.FALSE.)!
      PARAMETER       (YINTST=.TRUE.)
! Define 'backslash', some UNIX-systems want \\, some only \ ...
! so take the following and it should be highly portable ...
      PARAMETER       (BS='\\')
      EXTERNAL        LENGTH,NWORDS,NXTFRU,NITEMS

! Initialise N and IERR:
      N=0
      IERR=0

! Invalid data type given: return error code '4'
      IF ((TYPE.NE.'I').AND.(TYPE.NE.'F').AND. &
     &    (TYPE.NE.'C').AND.(TYPE.NE.'S').AND.(TYPE.NE.'L')) THEN
         IERR=4
         RETURN
      ENDIF

! Search some free unit where to open the input file FNAM if necessary:
      IF ((IU.LT.0).OR.(IU.GT.99)) IU=NXTFRU()
! Hmmm ... . No free unit found -- you should not use so much files!
      IF (IU.LT.0) THEN
         IERR=1
         RETURN
      ENDIF

! Try to open the data file ...
      IF (LOPEN) THEN
         WORK=FNAM
         CALL STRIP(WORK,L,'A')
         IF (L.EQ.0) GOTO 10
         OPEN(IU,FILE=WORK(1:L),STATUS='OLD',ERR=10)
      ELSE
         REWIND IU
      ENDIF
      GOTO 20
! Hmmm ... . No success? Maybe you should first create the file ...
   10 IERR=2
      IF (LOPEN) CLOSE(IU,ERR=11)
   11 CONTINUE
      RETURN
   20 CONTINUE
! Set some flag telling us whether we have found what we searched:
      FOUND=.FALSE.
! ... and set some flag telling us something about continuation lines:
      CONT=.FALSE.

! First we must somehow get information on the 'keyword' stored in WHAT:
! how many (non-blank) characters has the keyword really after stripping
! all leading and trailing blanks? Try to find out ...
      KEY=WHAT
      CALL STRIP(KEY,LKEY,'A')
      CALL UPPER(KEY)

! Now read the file line by line until we find the correct 'keyword'
! (or not ...) -- label 30 actually is an entry to an implicit loop!
   30 CONTINUE
! Read until the end of the file or until some error occurs ...
      READ(IU,'(A)',ERR=10000,END=10000) BUFLIN

! Continuation line has to be appended to current string!
      IF (CONT) THEN
         WORK=WORK(1:LMAX)//BUFLIN
         BUFLIN=WORK
      ENDIF
! Where is the 'end of the input string' ...:
      LMAX=LENGTH(BUFLIN)
! Continuation line(s) to be expected?
      INOSL=0
      DO 40 I=LMAX,1,-1
         IF (BUFLIN(I:I).NE.BS) THEN
            INOSL=I
            GOTO 50
         ENDIF
   40 CONTINUE
   50 CONT=(MOD(LMAX-INOSL+2,2).EQ.1)
! Read next line until this entry is complete ...:
      IF (CONT) THEN
         LMAX=LMAX-1
         WORK=BUFLIN(1:LMAX)
         GOTO 30
      ENDIF

! Forget about leading blanks ...:
      CALL STRIP(BUFLIN,LMAX,'L')
! The real end is given by character CCHAR: all behind it is a comment!
      IF (BUFLIN(1:1).EQ.CCHAR) THEN
         BUFLIN=' '
         LMAX=0
         GOTO 70
      ENDIF
      DO 60 I=2,LMAX
         IF ((BUFLIN(I:I).EQ.CCHAR).AND.(BUFLIN(I-1:I-1).NE.BS)) THEN
            BUFLIN=BUFLIN(1:I-1)
            LMAX=I-1
            GOTO 70
         ENDIF
   60 CONTINUE
   70 CONTINUE
! Following means: only blanks after removing comment, nothing is left:
      IF (LMAX.EQ.0) GOTO 30

! Do we have more than one definition on one line? ---> parse!
      IF (LMAX.EQ.1) THEN
         IF (BUFLIN(1:1).EQ.SCHAR) THEN
            WORK=' '
            BUFLIN=' '
            LMAX=0
            GOTO 90
         ENDIF
      ENDIF
      IF (BUFLIN(1:1).EQ.SCHAR) THEN
         WORK=' '
         BUFLIN=BUFLIN(2:LMAX)
         LMAX=LMAX-1
         GOTO 90
      ENDIF
      DO 80 I=2,LMAX
         IF ((BUFLIN(I:I).EQ.SCHAR).AND.(BUFLIN(I-1:I-1).NE.BS)) THEN
            WORK=BUFLIN(1:I-1)
            IF (I.LT.LMAX) BUFLIN=BUFLIN(I+1:LMAX)
            IF (I.EQ.LMAX) BUFLIN=' '
            LMAX=LMAX-I
            GOTO 90
         ENDIF
   80 CONTINUE
! Was apparently the last definition (nomore separators found ...):
      WORK=BUFLIN
      BUFLIN=' '
      LMAX=0
   90 WPARSE=WORK
! Okay! Is the 'keyword' there?
      CALL STRIP(WPARSE,L,'A')
! Well, some separator at the beginning of the line!? Try the next ...
      IF (L.EQ.0) GOTO 70
      CALL UPPER(WPARSE)
      IKEY=INDEX(WPARSE(1:L),KEY(1:LKEY))
! The keyword must be the first what appears in the definition string!!
! After stripping  all  blanks it must start at the first character ---
! otherwise we have found something else (but not meaningful!!!!!!!!!!):
      IF (IKEY.EQ.1) THEN
! Maybe we got it! But only if a parsing character follows immediately!
         IF ((1+LKEY).GT.L) THEN
! Impossible: PCHAR would be beyond the end of the string if this holds!
! Try the next definition ...
            GOTO 70
         ENDIF
         IF (WPARSE(1+LKEY:1+LKEY).NE.PCHAR) THEN
! No! We have found the word (substring ...) but no assignment character
! immediately after it --- must be some other entry or 'open comment'
! containing the keyword just as a substring; try next definition ...
            GOTO 70
         ENDIF
! Okay! Got it!!!!
         BUFLIN=WORK
         LMAX=LENGTH(BUFLIN)
         GOTO 100
      ELSE
! Try the next definition ...
         GOTO 70
      ENDIF
  100 FOUND=.TRUE.

! Parse off the data list now!
      CALL PARSE(BUFLIN,WPARSE,WORK,PCHAR,0)
      BUFLIN=WORK
      L=MAX(LENGTH(WPARSE),1)
! Hmmm ... . Have not yet found the correct parsing character, next try!
      IF (WPARSE(L:L).EQ.BS) GOTO 100

! And here we have it: get the data from the list!
      IF (TYPE.EQ.'S') THEN
! String -- that is tooooooooo easy ... :
         WORK=BUFLIN
! 'translate' all 'special characters':
         N=0
! the 'parsing character' ...
         HIDDEN=BS//PCHAR
         CALL REPLAC(WORK,HIDDEN,PCHAR,N,0)
! the 'comment character' ...
         HIDDEN=BS//CCHAR
         CALL REPLAC(WORK,HIDDEN,CCHAR,N,0)
! the 'separation character' ...
         HIDDEN=BS//SCHAR
         CALL REPLAC(WORK,HIDDEN,SCHAR,N,0)
! and the 'backslash' ...
         HIDDEN=BS//BS
         CALL REPLAC(WORK,HIDDEN,BS,N,0)
! You do not want more than NMAX characters ??
         NREAD=N
         IF (NMAX.GT.0) NREAD=MIN(N,NMAX)
         STR=WORK(1:NREAD)
      ELSE
! List of numbers/logicals (checking also types) --> how much valid??
         N=NITEMS(BUFLIN,WORK,.FALSE.,TYPE)
! There seems to be a fatal format error?
         IF (N.LT.1) GOTO 110
! Well, maybe we want/cant never read more than NMAX data ...
         NREAD=N
         IF (NMAX.GT.0) NREAD=MIN(N,NMAX)
! Now really read the data (no responsibilty for core dumps here due to
! insufficient dimensioning of the arrays 'coming from above' ... !!):
         IF (YINTST) THEN
! Well, here is the 'critical stuff' with YINTST ... . If set .FALSE. the
! compiler should hopefully note "code unreachable" and continue without
! complaining about the invalid '*'-format ... . If you still run into
! trouble because your compiler tries to translate the following section
! although it can never be executed then comment out the following block
! in order to get the code working ... (shit FORTRAN77 ... !).
            IF (TYPE.EQ.'I') THEN
               READ(BUFLIN,*,ERR=110,END=110) (INTRES(I),I=1,NREAD)
            ELSE IF (TYPE.EQ.'F') THEN
               READ(BUFLIN,*,ERR=110,END=110) (FLTRES(I),I=1,NREAD)
            ELSE IF (TYPE.EQ.'C') THEN
               READ(BUFLIN,*,ERR=110,END=110) (CMPRES(I),I=1,NREAD)
            ELSE IF (TYPE.EQ.'L') THEN
               READ(BUFLIN,*,ERR=110,END=110) (LOGRES(I),I=1,NREAD)
            ENDIF
         ELSE
            IUSCR=NXTFRU()
            IF (IUSCR.LT.0) THEN
               IERR=6
               N=0
               RETURN
            ENDIF
            OPEN(IUSCR,STATUS='SCRATCH',ERR=130)
            WRITE(IUSCR,'(A)') BUFLIN
            REWIND IUSCR
            IF (TYPE.EQ.'I') THEN
               READ(IUSCR,*,ERR=110,END=110) (INTRES(I),I=1,NREAD)
            ELSE IF (TYPE.EQ.'F') THEN
               READ(IUSCR,*,ERR=110,END=110) (FLTRES(I),I=1,NREAD)
            ELSE IF (TYPE.EQ.'C') THEN
               READ(IUSCR,*,ERR=110,END=110) (CMPRES(I),I=1,NREAD)
            ELSE IF (TYPE.EQ.'L') THEN
               READ(IUSCR,*,ERR=110,END=110) (LOGRES(I),I=1,NREAD)
            ENDIF
            CLOSE(IUSCR,ERR=130)
         ENDIF
         GOTO 120
! Hmmm ... . Something was not okay with the format ...
  110    IERR=5
         N=0
         GOTO 10000
  120    CONTINUE
         GOTO 140
  130    IERR=6
         N=0
         CLOSE(IUSCR,ERR=10000)
         GOTO 10000
  140    CONTINUE
      ENDIF

! Yepeeh! That is the end folks -- successful or not (who cares ...)!
10000 CONTINUE
! Close the input file ...
      IF (LOPEN) CLOSE(IU,ERR=10001)
! Hmmm ... . No success? Maybe you should add something in your data
! file, or you should check the format ... (??)  --  LOOSER!!
      IF (.NOT.FOUND) IERR=3
10001 CONTINUE
! ... and bye bye my honey -- would be nice to meet you again ...
      RETURN

      END
