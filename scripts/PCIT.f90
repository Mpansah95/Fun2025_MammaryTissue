!###############################################################################
!                                                                              #
!         CCCCCC  SSSSSS  IIIII  RRRRRR   OOOO          L      IIIII           #
!         C       S         I    R    R  O    O         L        I             #
!         C       SSSSSS    I    RRRRRR  O    O  =====  L        I             #
!         C            S    I    R RR    O    O         L        I             #
!         CCCCCC  SSSSSS  IIIII  R   RR   OOOO          LLLLL  IIIII           #
!                                                                              #
!         Queensland Bioscience Precint                                        #
!         306 Carmody Road                                                     #
!         St Lucia, QLD 4067                                                   #
!         AUSTRALIA                                                            #
!         Ph. (07) 3214 2392                                                   #
!         Fx. (07) 3214 2900                                                   #
!                                                                              #
!         Objective:      Perform PCIT algorithm to build gene networks        #
!         Authors:        Antonio Reverter (Tony.Reverter-Gomez@csiro.au)      #
!                         Eva Chan (Eva.Chan@csiro.au)                         #
!         Last Revision:  January 2011                                         #
!                                                                              #
!###############################################################################
PROGRAM pcit_simulateddata

IMPLICIT NONE

INTEGER  :: ngenes  ,&     !=4120, &  ! Number of genes in the network
                      nconditions    !=40     ! Experimenal conditions = dimension
                                             ! of gene expression vectors

TYPE gene_record
   CHARACTER(LEN=40) :: id
   !double precision :: expr(nconditions)
   double precision,allocatable :: expr(:)
END TYPE gene_record
!TYPE(gene_record) :: gene(ngenes)
TYPE(gene_record),allocatable :: gene(:)

!double precision :: gene_corr(ngenes,ngenes),  &
double precision ,allocatable:: gene_corr(:,:),  &
        gene_pcorr(:,:)
        !gene_pcorr(ngenes,ngenes), &
double precision::  val,ave,mini,maxi,corr,    &
        rxy, rxz, ryz,             &  ! Direct correlations
        rxy_g_z, rxz_g_y, ryz_g_x, &  ! Partial correlations
        toler                         ! Tolerance threshold

INTEGER :: i,j,k,io,maxlengthid
character(len=100):: infile, outfile, line,fmat

write(*,*) 'name of file with expressions for each gene across all samples/conditions'
read(*,*)infile
write(*,*)infile
outfile=trim(adjustl(infile))//'_out.txt'
write(*,*) 'the output file is ',trim(outfile)

OPEN(20,file=infile,status='old',form='formatted')
print *,1
OPEN(50,file=outfile,status='replace',form='formatted')
print *,2
!#################################################################


ngenes=0
do
  read(20,'(a1)',iostat=io)line(1:1)
  if(io/=0) exit
  ngenes=ngenes+1
enddo
write(*,*)'number of genes= ',ngenes
rewind(20)
call number_of_fields(20,nconditions)
nconditions=nconditions-1
write(*,*) 'number of conditions: ',nconditions

allocate(gene_corr(ngenes,ngenes),gene_pcorr(ngenes,ngenes), gene(ngenes))
gene_corr = 1.d-4
gene_pcorr = 1.d-4
do i=1,ngenes
    allocate(gene(i)%expr(nconditions))
enddo
maxlengthid=0
DO i = 1,ngenes
   !READ(20,'(a40,40f10.5)')gene(i)%id,gene(i)%expr(1:nconditions)
   READ(20,*)gene(i)%id,gene(i)%expr(1:nconditions)
   maxlengthid=maxval((/maxlengthid,len_trim(adjustl(gene(i)%id))/))
ENDDO
print *,'max length of gene name: ',maxlengthid

!###################################################
! Process Genes
!###################################################
DO i = 1, ngenes-1
   DO j = i+1, ngenes
      gene_corr(i,j) = corr(gene(i)%expr,gene(j)%expr,nconditions)
      gene_corr(j,i) = gene_corr(i,j)
      if(gene_corr(i,j)==0) then
          print *,i,j,gene_corr(i,j)
      endif
   ENDDO
ENDDO
!DO i = 1, ngenes
!   ave = SUM(gene_corr(i,:))/(ngenes-1)
!   mini = MINVAL(gene_corr(i,:))
!   maxi = MAXVAL(gene_corr(i,:))
!  WRITE(40,'(i7,i7,3f9.3)')i,ngenes-1,ave,mini,maxi
!ENDDO

!########################################################
! Get all combinations of genes taking 3 at a time
! For each trio, compute the 3 potential comparisons
! to stablish which correlations is to be cancelled out
!########################################################
gene_pcorr = gene_corr
DO i = 1, ngenes-2
   if(mod(i,10)==0)write(*,*)'Trios for gene',i
   DO j = i+1, ngenes-1
      DO k = j+1, ngenes
         rxy = gene_corr(i,j)
         rxz = gene_corr(i,k)
         ryz = gene_corr(j,k)

         rxy_g_z = (rxy - rxz*ryz) / sqrt((1-rxz**2) * (1-ryz**2))
         rxz_g_y = (rxz - rxy*ryz) / sqrt((1-rxy**2) * (1-ryz**2))
         ryz_g_x = (ryz - rxy*rxz) / sqrt((1-rxy**2) * (1-rxz**2))
         toler = abs(rxy_g_z/rxy)
         toler = toler + abs(rxz_g_y/rxz)
         toler = toler + abs(ryz_g_x/ryz)
         toler = 0.33333 * toler
         !WRITE(60,*)i,j,k,toler
         ! NB: This "Tolerance.out" output file could be huge
         !     with many genes. Consider not to produce it unless
         !     there is a major interest in looking at the tolerance
         !     threshold levels
         IF( abs(rxy) < abs(rxz*toler) .AND. abs(rxy) < abs(ryz*toler) )THEN
            gene_pcorr(i,j) = 0.
         ENDIF
         IF( abs(rxz) < abs(rxy*toler) .AND. abs(rxz) < abs(ryz*toler) )THEN
            gene_pcorr(i,k) = 0.
         ENDIF
         IF( abs(ryz) < abs(rxy*toler) .AND. abs(ryz) < abs(rxz*toler) )THEN
            gene_pcorr(j,k) = 0.
         ENDIF
      ENDDO
   ENDDO
ENDDO

! create format for output
write(fmat,'(a,i0,a)') '(2a',maxlengthid+1,',1x,f10.5,1x,f10.5)'
DO i = 1, ngenes-1
   DO j = i+1, ngenes
      gene_pcorr(j,i) = gene_pcorr(i,j)
      WRITE(50,fmat)gene(i)%id,gene(j)%id, &
                               gene_corr(i,j),gene_pcorr(i,j)
   ENDDO
ENDDO


END PROGRAM pcit_simulateddata

!#########################################################################
!#########################################################################
FUNCTION corr(x,y,n) RESULT(f_val)
IMPLICIT NONE

INTEGER, INTENT(IN) :: n
double precision, INTENT(IN) :: x(n),y(n)
double precision :: f_val

double precision :: s1, ss1, s2, ss2, ss12, var1, var2, num, den
INTEGER :: i,j

s1=0.0; ss1=0.0; s2=0.0; ss2=0.0; ss12=0.0

DO i = 1, n
   s1 = s1 + x(i)
   ss1 = ss1 + x(i)*x(i)
   s2 = s2 + y(i)
   ss2 = ss2 + y(i)*y(i)
   ss12 = ss12 + x(i)*y(i)
ENDDO
var1 = (ss1 - (s1*s1)/n) / (n-1)
var2 = (ss2 - (s2*s2)/n) / (n-1)
num = (ss12 - (s1*s2)/n) / (n-1)
den = SQRT(var1*var2)

f_val = num / den
if(f_val==0) then 
    print *,num,den,num/den
endif
END FUNCTION corr


subroutine number_of_fields(unit,n)
    ! guess number of fields in the 1st line of one file
 integer unit,n,i
 character(100000):: line
 character(1)::last
 last=' '
 n=0
 read(unit,'(a)')line
 do i=1,len_trim(line)
    if (line(i:i)/=' ' .and. last==' ') then
        n=n+1
    endif
    last=line(i:i)
 enddo
 !write (*,*) n
 rewind(unit)
end subroutine 

