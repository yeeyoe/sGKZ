#import <Foundation/Foundation.h>

NS_ASSUME_NONNULL_BEGIN

@interface GKZGuiResult : NSObject
@property(nonatomic, readonly) BOOL certified;
@property(nonatomic, readonly) BOOL cancelled;
@property(nonatomic, readonly) NSString *message;
@property(nonatomic, readonly) NSArray<NSString *> *sigmaExact;
@property(nonatomic, readonly) NSArray<NSString *> *sigmaVeeExact;
@property(nonatomic, readonly) NSArray<NSArray<NSNumber *> *> *triangulationFaces;
@property(nonatomic, readonly) NSArray<NSArray<NSValue *> *> *subdivisionCells;
@property(nonatomic, readonly) NSInteger iterations;
@property(nonatomic, readonly) NSUInteger activeSize;
@end

NS_ASSUME_NONNULL_END
