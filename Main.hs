import Numeric.LinearAlgebra
import Numeric.LinearAlgebra.HMatrix
import Prelude hiding (pi)

type X = [V]
type V = Vector Double
type M = Matrix Double

data Theta = Theta {mu::[V], sigma::[M], pi::[Double]}
    deriving Show

numCluster = 4
dim = 2
numDataset = 3
numIterate = 10

vec = vector [1,2]

initTheta :: Theta
initTheta = Theta [vector [0.1,0.1],vector [0.1,0.8],vector [0.8,0.1]] [i,i,i] [0.1,0.2,0.7]
    where i = (2><2) [1.0,0,0,0,1.0,0,0,0,1.0]
initX :: X
initX = [vector [0.5,0.5],vector [0.2,0.5],vector [0.3,0.3]]

-- for given dataset and initial parameter return list of cluster for each data
em :: X -> Theta -> [[Double]]
em x t = [[gamma x t n k|k<-[0..numCluster-1]]|n<-[0..numDataset-1]]
emParameter = (iterate (mStep X) t) !! numIterate

eStep :: X -> Theta -> M
eStep x t = (numDataset><numCluster) (concat [[gamma x t n k|k <- [0..numCluster-1]]|n <- [0..numDataset-1]])
-- gamma _ _ _ _ = 1
gamma x t n k = ((pi t !! k) * gauss (x !! n) (mu t !! k) (sigma t !! k)) / sum [((pi t !! j) * gauss (x !! n) (mu t !! j) (sigma t !! j))|j<-[0..numCluster-1]]

mStep :: X -> Theta -> Theta
mStep x t = Theta muNew sigmaNew piNew
    where n :: Int -> Double
          n k = sum [gamma x t i k|i<-[0..numDataset-1]]
          muNew :: [V]
          muNew = [sum [scalar ((gamma x t i k)/(n k))*(x !! i)|i<-[0..numDataset]]|k<-[0..numCluster]]
          sigmaNew :: [M]
          sigmaNew = [sum [scalar (gamma x t i k/n k)*mul (a i k) (b i k)|i<-[0..numDataset]]|k<-[0..numCluster]]
          a i k = matrix dim (toList ((x !! i) - (muNew !! k)))
          b i k = matrix 1 (toList ((x !! i) - (muNew !! k))) 
          piNew :: [Double]
          piNew = [n k/fromIntegral numDataset|k<-[0..numCluster]]

lnLikelihood :: X -> Theta -> Double
lnLikelihood x t = sum [ log (sum [ (pi t !! k) * gauss (x !! n) (mu t !! k) (sigma t !! k) | k<-[0..numCluster-1] ]) | n<-[0..numDataset-1] ]

-- gauss x mu sigma 
gauss :: V -> V -> M -> Double
gauss x m s = exp (((asColumn y) <> t <> (asRow y)) ! 1 ! 1) / ((sqrt 2*3.14)^dim*abs (det s))
    where y = x - m
          t = inv s

main = putStrLn "Hello"

