{-# LANGUAGE ScopedTypeVariables #-}

module Dna.Graph (
readConnection,readAllConnection,allPossibleNode,graphFromConnection,Graph,euclidianCycle
)
where

import Data.List
import Data.Functor
import Data.STRef
import Control.Monad
import Data.Maybe
import Control.Applicative
import Numeric.LinearAlgebra (rows, (!), cols, toLists, (?), toList)
import Numeric.LinearAlgebra.Data (Matrix, I)
import Numeric.LinearAlgebra.Devel
import Control.Monad.ST

import qualified Data.ByteString.Lazy.Char8 as C
import qualified Data.Set as S

data Ord a => Graph a = GraphInternal {
    node :: S.Set a,
    connection :: Matrix I
}

instance (Ord a, Show a) => Show (Graph a) where
    show graph = intercalate "\n"
            [link | row <- [0 .. rows (connection graph) - 1],
                    any (isLink row) [0 .. cols (connection graph) - 1],
                    link <- [show (S.elemAt row (node graph)) ++ " -> " ++
                             intercalate "," (concatMap (\col -> replicate (fromIntegral $ linkCount row col) $ show  $ S.elemAt col (node graph))
                                                (filter (isLink row) [0 .. cols (connection graph) - 1]))]]
                where
                    linkCount row col = connection graph ! row ! col
                    isLink row col =   linkCount row col >= 1

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
    where
        nodes = S.fromList $ allPossibleNode connections
        findIndex = flip S.findIndex nodes
        nodesLen = S.size nodes

euclidianCycle :: Ord a => Graph a -> Maybe [(a,a)]
euclidianCycle (GraphInternal nodes m) = runST $ do
    traverseM <- thawMatrix m
    path <- findEuclidianCycleAt nodes traverseM (Just 0)
    findNext path traverseM
        where findNext path m = do
                newPath <- findNextEuclidianPath nodes m path
                if isJust newPath
                    then findNext newPath m
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


euclidianNode :: Int -> STMatrix s I -> ST s Bool
euclidianNode n m = do
    colm <- extractMatrix m AllRows (Col n)
    rowm <- extractMatrix m (Row n) AllCols
    return $ (sum . concat . toLists) colm == (sum . concat . toLists) rowm
