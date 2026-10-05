#import <Foundation/Foundation.h>
#import "GKZGuiResult.h"

NS_ASSUME_NONNULL_BEGIN

typedef void (^GKZProgressBlock)(NSInteger iterations, NSUInteger activeSize,
                                 NSString *status);

@interface GKZCancellationToken : NSObject
- (void)cancel;
@property(nonatomic, readonly, getter=isCancelled) BOOL cancelled;
@end

@interface GKZBridge : NSObject
+ (GKZGuiResult *)solvePoints:(NSArray<NSValue *> *)points
                     progress:(nullable GKZProgressBlock)progress
                        token:(GKZCancellationToken *)token;
@end

NS_ASSUME_NONNULL_END
