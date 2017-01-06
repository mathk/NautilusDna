{-# LANGUAGE ScopedTypeVariables, RankNTypes #-}

module Dna.Graph (
readConnection,readAllConnection,
allPossibleNode,graphFromConnection,
Graph,euclidianCycle,
unbalanceNodes,
canBeBalance, balance, initOverlaGraph
)
where

import Data.List
import Data.Tuple
import Data.Functor
import Data.STRef
import Control.Monad
import Data.Maybe
import Control.Applicative
import Numeric.LinearAlgebra (rows, (!), cols, toLists, (?), (¿), toList)
import Numeric.LinearAlgebra.Data (Matrix, I)
import Numeric.LinearAlgebra.Devel
import Control.Monad.ST

import qualified Data.ByteString.Lazy.Char8 as C
import qualified Data.Set as S

data Graph a = GraphInternal {
    node :: S.Set a,
    connection :: Matrix I,
    addCon :: Maybe (a,a)
}

instance (Ord a, Show a) => Show (Graph a) where
    show graph = fromMaybe matrixShow (do
                    (from,to) <- addCon graph
                    return $ matrixShow ++ "\n With Link: " ++ show from  ++ "->" ++ show to)
        where
            linkCount row col = connection graph ! row ! col
            isLink row col =   linkCount row col >= 1
            matrixShow = intercalate "\n"
                [link | row <- [0 .. rows (connection graph) - 1],
                        any (isLink row) [0 .. cols (connection graph) - 1],
                        link <- [show (S.elemAt row (node graph)) ++ " -> " ++
                                 intercalate "," (concatMap (\col -> replicate (fromIntegral $ linkCount row col) $ show  $ S.elemAt col (node graph))
                                                    (filter (isLink row) [0 .. cols (connection graph) - 1]))]]

initOverlaGraph :: Ord a =>  [a] -> (forall s. STMatrix s I -> (a -> Int) -> ST s ()) -> Graph a
initOverlaGraph nodes f = GraphInternal setNode (runSTMatrix $ do
        m <- newMatrix 0 lenNodes lenNodes
        f m (`S.findIndex` setNode)
        return m
    ) Nothing
    where
        setNode = S.fromList nodes
        lenNodes  = length setNode


readConnection :: (Read a) => String -> (a, [a])
readConnection info = case lex info of
    [(node, left)] -> (read node, readListOfNode (snd $ head $ lex left))

readListOfNode :: (Read a) => String -> [a]
readListOfNode nodes = case lex nodes of
    [("", "")] -> []
    [(node, left)] -> read node : readListOfNode (snd $ head $ lex left)

readAllConnection :: (Read a) => String -> [(a, [a])]
readAllConnection = map readConnection . lines

allPossibleNode :: [(a,[a])] -> [a]
allPossibleNode = concatMap $ uncurry (:)

graphFromConnection ::Ord a =>  [(a,[a])] -> Graph a
graphFromConnection connections = GraphInternal
        nodes
        (runSTMatrix $ do
            m <- newMatrix 0 nodesLen nodesLen
            forM_ connections $ \ (from,connectee) ->
                forM_ connectee $ \ to ->
                    modifyMatrix m (findIndex from) (findIndex to) (1+)
            return m)
        Nothing
    where
        nodes = S.fromList $ allPossibleNode connections
        findIndex = flip S.findIndex nodes
        nodesLen = S.size nodes

euclidianCycle :: Ord a => Graph a -> Maybe [(a,a)]
euclidianCycle g = do
    (GraphInternal nodes m newCon) <- balance g
    conn <- runST $ do
        traverseM <- thawMatrix m
        path <- findEuclidianCycleAt nodes traverseM (Just 0)
        findNext nodes path traverseM
    return $ maybe conn (\(from,to) -> init $ uncurry (++) (swap (break ( (to == ).fst ) conn))) newCon
  where findNext nodes path m = do
            newPath <- findNextEuclidianPath nodes m path
            if isJust newPath
                then findNext nodes newPath m
                else return path

findNextEuclidianPath :: Ord a => S.Set a -> STMatrix s I -> Maybe [(a,a)] -> ST s (Maybe [(a,a)])
findNextEuclidianPath nodes mat path = do
    inMat <- freezeMatrix mat
    let newPath = do
            foundPath <- path
            let (left,right) = break (\(from,to) -> sum (map ti $ toList (inMat ! S.findIndex from nodes)) > 0) foundPath
            if null right
                then Nothing
                else Just (right ++ left)
    let maybeNext = fmap ((`S.findIndex` nodes) . fst . head) newPath
    maybeNextPath <- findEuclidianCycleAt nodes mat maybeNext
    return $ liftM2 (++) maybeNextPath newPath


findEuclidianCycleAt :: Ord a => S.Set a -> STMatrix s I -> Maybe Int -> ST s (Maybe [(a,a)])
findEuclidianCycleAt nodes mat startNode = if isJust startNode
                                                then findNextNod startNode
                                                else return Nothing
    where
        findNextNod n =  do
            next <- findNextEuclidianAt nodes mat (fromJust n)
            let cycleEnd = liftM2 (==) next startNode
            case (cycleEnd, next) of
                (Just True, Just to) -> return (Just [(S.elemAt (fromJust n) nodes,S.elemAt to nodes)])
                (Just False, to) -> do
                    list <- findNextNod to
                    return $ liftM2 (:) (Just (S.elemAt (fromJust n) nodes,S.elemAt (fromJust to) nodes)) list
                _ -> return Nothing

findNextEuclidianAt :: Ord a => S.Set a -> STMatrix s I -> Int -> ST s (Maybe Int)
findNextEuclidianAt nodes mat n = do
    con <- mapM (readMatrix mat n) [0..(S.size nodes - 1)]
    let next = find ((>0).snd) (zip [0..] con)
    when (isJust next) $ modifyMatrix mat n (fst . fromJust $ next) (subtract 1)
    return $ fst <$> next

mapFst :: (a -> b) -> [(a,c)] -> [(b,c)]
mapFst f = map (uncurry $ (,) . f)

canBeBalance :: Ord a => Graph a -> Bool
canBeBalance g = case unbal of
    [(_,-1), (_,1)] -> True
    [(_,1), (_,-1)] -> True
    [] -> True
    _ -> False
  where unbal = unbalanceNodes g

balance :: Ord a => Graph a -> Maybe (Graph a)
balance g@(GraphInternal node connection addCon) = case unbal of
    [(to,-1), (from,1)] -> return (GraphInternal node (fst (mutable (correctingMatrix (S.findIndex from node) (S.findIndex to node)) connection)) (Just (from,to)))
    [(from,1), (to,-1)] -> return (GraphInternal node (fst (mutable (correctingMatrix (S.findIndex from node) (S.findIndex to node)) connection)) (Just (from,to)))
    [] -> return g
    _ -> Nothing
  where
    unbal = unbalanceNodes g
    correctingMatrix f t _ m = modifyMatrix m f t (+1)

unbalanceNodes :: Ord a => Graph a -> [(a,Int)]
unbalanceNodes (GraphInternal nodes m _) =filter ((0 /= ) . snd) (mapFst (`S.elemAt` nodes) (map (\idx -> (idx,eulerianNode idx m)) [0..(S.size nodes - 1)]))

eulerianNode :: Int -> Matrix I -> Int
eulerianNode n m = ti (sum (concat colConnection)) - ti (sum (concat rowConnection))
    where
        colConnection = toLists $ m ¿ [n]
        rowConnection = toLists $ m ? [n]

eulerianNodeST :: Int -> STMatrix s I -> ST s Bool
eulerianNodeST n m = do
    colm <- extractMatrix m AllRows (Col n)
    rowm <- extractMatrix m (Row n) AllCols
    return $ (sum . concat . toLists) colm == (sum . concat . toLists) rowm
