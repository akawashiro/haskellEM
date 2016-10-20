import Numeric.LinearAlgebra
import Numeric.LinearAlgebra.HMatrix
import Prelude hiding (pi)
import Text.Printf

type X = [V]
type V = Vector Double
type M = Matrix Double

data Theta = Theta {mu::[V], sigma::[M], pi::[Double]}
    deriving Show

numCluster = 4
dim = 2
numDataset = length initX
numIterate = 10

vec1 = vector [1,2]
vec2 = vector [2,1]
mat1 :: M
mat1 = (2><2) [1.0,0,0,1.0]

initTheta :: Theta
initTheta = Theta [vector [0.2,0.2],vector [0.2,0.8],vector [0.8,0.2],vector [0.8,0.8]] [i,i,i,i] [0.25,0.25,0.25,0.25]
    where i = (2><2) [1,0,0,1]

initX :: X
initX = concat [cluster1,cluster2,cluster3,cluster4,cluster5]
    where
        clusterFromList xs ys = concat [[[x,y]|x<-xs]|y<-ys]
        cluster1 = map vector $ clusterFromList [0.1,0.2,0.3] [0.1,0.2,0.3]
        cluster2 = map vector $ clusterFromList [0.1,0.2,0.3] [0.7,0.8,0.9]
        cluster3 = map vector $ clusterFromList [0.7,0.8,0.9] [0.1,0.2,0.3]
        cluster4 = map vector $ clusterFromList [0.7,0.8,0.9] [0.7,0.8,0.9]
        cluster5 = map vector $ clusterFromList [0.4,0.5,0.6] [0.4,0.5,0.6]

testEM = em initX initTheta

-- for given dataset and initial parameter return list of cluster for each data
em :: X -> Theta -> [[Double]]
em x t = [[gamma x (emParameter x t) n k|k<-[0..numCluster-1]]|n<-[0..numDataset-1]]
emParameter x t = (iterate (mStep x) t) !! numIterate

eStep :: X -> Theta -> M
eStep x t = (numDataset><numCluster) (concat [[gamma x t n k|k <- [0..numCluster-1]]|n <- [0..numDataset-1]])
gamma x t n k = ((pi t !! k) * gauss (x !! n) (mu t !! k) (sigma t !! k)) / sum [((pi t !! j) * gauss (x !! n) (mu t !! j) (sigma t !! j))|j<-[0..numCluster-1]]

mStep :: X -> Theta -> Theta
mStep x t = Theta muNew sigmaNew piNew
    where n :: Int -> Double
          n k = max 0.0000001 (sum [gamma x t i k|i<-[0..numDataset-1]])
          muNew :: [V]
          muNew = [sum [scalar ((gamma x t i k)/(n k))*(x !! i)|i<-[0..numDataset-1]]|k<-[0..numCluster-1]]
          sigmaNew :: [M]
          sigmaNew = [sum [scalar (gamma x t i k/n k)*c i k|i<-[0..numDataset-1]]|k<-[0..numCluster-1]]
          a i k = asColumn ((x !! i) - (muNew !! k))
          b i k = asRow ((x !! i) - (muNew !! k))
          c i k = mul (a i k) (b i k)
          piNew :: [Double]
          piNew = [n k/fromIntegral numDataset|k<-[0..numCluster-1]]

lnLikelihood :: X -> Theta -> Double
lnLikelihood x t = sum [ log (sum [ (pi t !! k) * gauss (x !! n) (mu t !! k) (sigma t !! k) | k<-[0..numCluster-1] ]) | n<-[0..numDataset-1] ]

-- gauss x mu sigma 
gauss :: V -> V -> M -> Double
gauss x m s = exp (-1*((asRow y) <> t <> (asColumn y)) ! 0 ! 0) / ((sqrt 2*3.14)^dim*abs (det s))
    where y = x - m
          t = inv s

main = do
    putStr "Data points are \n"
    putStr (show initX)
    putStr "\n\n"
    putStr "Initail parameters are \n"
    putStr (show initTheta)
    putStr "\n\n"
    putStr "After EM algorithm each data points are classified to cluster following probabilities\n"
    putStr "[probability of cluster 0, probability of cluster1, probability of cluster2, probability of cluster3]\n\n"
    putStrLn ("[" ++ (showResult testEM))

showResult :: [[Double]] -> String
showResult [] = "]"
showResult (r:rs) = "[" ++ showResult' r ++ "], \n" ++ showResult rs
showResult' [] = ""
showResult' [x] = printf "%.2f" x
showResult' (x:xs) = printf "%.2f," x ++ showResult' xs

